function add_awgn(X::AbstractMatrix, snr_db::Real)
    P = mean(abs2, X)
    X .+ sqrt(P / 10^(snr_db / 10)) .* randn(eltype(X), size(X))
end

"""Double-exponential calcium transient kernel (GCaMP-like)."""
function make_soma_activity(K, nt, n_spikes;
                             tau_rise=3f0, tau_decay=60f0, kernel_mult=6, amp=4f0, seed=nothing)
    seed === nothing ? Random.seed!() : Random.seed!(seed)
    kernel_len = round(Int, kernel_mult * tau_decay)
    ker = Float32.([exp(-t / tau_decay) - exp(-t / tau_rise) for t in 0:kernel_len-1])
    ker ./= maximum(ker)
    soma = fill(1f0, K, nt)
    for k in 1:K
        times = sort(randperm(nt - 80)[1:n_spikes] .+ 40)
        for t0 in times
            len = min(length(ker), nt - t0 + 1)
            soma[k, t0:t0+len-1] .+= amp .* ker[1:len]
        end
    end
    soma
end

"""Soma activity from NAOMiSim's own spike-generation + calcium-dynamics
pipeline (Poisson/burst spiking → Ca_DE binding kinetics), with baseline
fluorescence individualized per neuron via NAOMiSim's built-in expression
variation (`SpikeOpts.min_mod`/`p_off`, see `expression_variation`).

Spike generation runs internally at 100 Hz regardless of `dt`
(`generate_time_traces`), so `rate` (events per internal 100 Hz frame) must be
calibrated against that internal frame count, not against `nt`/`dt` directly
— hence the extra `/100`. `smod_flag="hawkes"` (NAOMiSim's default) is
avoided: with `SpikeOpts`'s default self-excitation (`selfact=1.2`, supercritical)
the Hawkes process saturates near its ceiling almost continuously regardless
of `rate`, washing out the transient pulses; plain Poisson/burst spiking
(`smod_flag="burst"`) stays close to baseline between bursts instead.

NAOMiSim has no concept of inhibitory neurons (all spikes raise fluorescence).
`inh_frac` designates a random fraction of neurons as inhibitory by mirroring
their trace about its own per-neuron baseline (`mod_vals[k]`) in log/ratio
space (`baseline^2 / x`), so spikes *decrease* fluorescence for those neurons
instead of increasing it. A linear (additive) mirror would drive most of the
trace negative — GCaMP excursions reach tens of times baseline — requiring
clamping that flattens the pulse shape; the ratio form stays positive and
keeps the dip shape smooth. Inhibitory neurons also get their resting
baseline raised by `inh_baseline_boost`, so the dips have headroom and don't
crowd against zero."""
function make_naomi_soma_activity(K, nt, dt, n_spikes;
                                   prot="GCaMP6", dyn_type="Ca_DE",
                                   rate_dist="gamma", alpha=1.0,
                                   smod_flag="burst", min_mod=[0.4, 2.53], p_off=0.0,
                                   inh_frac=0.0, inh_baseline_boost=2.0,
                                   seed=nothing)
    seed === nothing ? Random.seed!() : Random.seed!(seed)
    rate = n_spikes / (nt * 100 * dt)   # events per internal 100 Hz frame
    spike_opts = SpikeOpts(; K, nt, dt, rate, rate_dist, alpha, dyn_type, prot,
                            min_mod, p_off, smod_flag, N_bg=0, dendflag=false, axonflag=false)
    S, _, mod_vals = generate_time_traces(spike_opts) # mod_vals are per-neuron baselines (expression levels)
    soma   = S.soma
    inh_idx = randperm(K)[1:round(Int, inh_frac * K)]
    for k in inh_idx
        soma[k, :] .= inh_baseline_boost .* mod_vals[k]^2 ./ soma[k, :]
        mod_vals[k] *= inh_baseline_boost
    end
    soma, mod_vals, inh_idx
end

"""Set nucleus voxels to nuc_val in gp_vals (makes dark nuclear hole)."""
function apply_nucleus_darkening!(vol_out, nuc_val=0f0)
    K = size(vol_out.locs, 1)
    for kk in 1:K
        gp_nuc_kk = vol_out.gp_nuc[kk]
        isnothing(gp_nuc_kk) && continue
        nuc_idxs = gp_nuc_kk[1]
        isempty(nuc_idxs) && continue
        nuc_set  = Set(nuc_idxs)
        soma_idxs = vol_out.gp_vals[kk][1]
        soma_vals = vol_out.gp_vals[kk][2]
        for ii in eachindex(soma_idxs)
            soma_idxs[ii] in nuc_set && (soma_vals[ii] = nuc_val)
        end
    end
end

"""Set all soma voxel fluorescence to soma_val (uniform cytoplasm)."""
function make_uniform_soma!(vol_out, soma_val=1f0, nuc_val=0f0)
    for kk in 1:size(vol_out.locs, 1)
        vol_out.gp_vals[kk][2] .= soma_val
    end
    apply_nucleus_darkening!(vol_out, nuc_val)
end

"""Add a white border around each panel image."""
function add_border(img::AbstractMatrix, bw::Int=1, val=maximum(img))
    h, w = size(img)
    out  = fill(Float64(val), h + 2bw, w + 2bw)
    out[bw+1:bw+h, bw+1:bw+w] .= img
    out
end

"""Save a 2D array as a grayscale PNG (without going through Makie)."""
function save_gray_png(path, A::AbstractMatrix, lo, hi; scale::Int=1)
    norm_A = clamp.((A .- lo) ./ (hi - lo), 0.0, 1.0)
    img = Gray.(N0f8.(norm_A))
    scale > 1 && (img = repeat(img, inner=(scale, scale)))
    save(path, img)
end

"""Save selected frames of a (rows, cols, time) array as an animated grayscale GIF
(without going through Makie)."""
function save_gray_gif(path, vol::AbstractArray{<:Real,3}, frames, lo, hi; scale::Int=1, fps::Int=8)
    norm_vol = clamp.((vol[:, :, frames] .- lo) ./ (hi - lo), 0.0, 1.0)
    img = Gray.(N0f8.(norm_vol))
    scale > 1 && (img = repeat(img, inner=(scale, scale, 1)))
    save(path, img; fps)
end

#== Prescan: generate and save ground-truth W, H ==#
function naomi_prescan(params, noise_params, tpm_params, imgsz, verbose)
    @unpack seed, N_neur, vol_sz, vol_depth, min_dist, avg_rad, nuc_rad, vres,
    nt, dt, prot, vasc_flag, psf_type, psf_NA, scan_buff, sfrac, n_spikes, spike_amp,
    inh_frac, inh_baseline_boost = params
    N1, N2 = imgsz

    # NAOMiSim has no seed parameter of its own (vol/dendrite generation use the
    # global RNG), so seed it here to make the whole volume reproducible.
    seed === nothing ? Random.seed!() : Random.seed!(seed)

    neur_params = NeurParams(;
        avg_rad     = avg_rad,
        nuc_rad     = nuc_rad,
        nuc_fluorsc = 0.0,       # nucleus has no fluorescence (GCaMP is cytoplasmic)
    )
    dend_params = DendParams(    # minimize dendrites to avoid grid artifacts
        dtParams  = [0., 2., 2., 1., 0.],
        atParams  = [0., 0., 0., 0., 0.],
        atParams2 = [0., 0., 0., 0., 0.],
        dweight   = 0.0,
    )
    vol_params  = VolumeParams(; vol_sz, vol_depth, N_neur, min_dist, vres,
                                    AD_density=0., verbose=1)
    psf_params  = PSFParams(; type=psf_type, NA=psf_NA)
    vasc_params = VascParams(; flag=vasc_flag)
    bg_params   = BgParams(; flag=false)
    axon_params = AxonParams(; flag=false)

    @info "Simulating neural volume ($(N_neur) neurons, $(vol_sz) μm)…"
    vol_out, vol_params = simulate_neural_volume(
        vol_params, neur_params, vasc_params, dend_params, bg_params, axon_params)[1:2]

    K = size(vol_out.locs, 1)
    @info "Done — $K neurons  z=$(round.(vol_out.locs[:,3]; digits=1)) μm"

    # Uniform cytoplasm fluorescence + dark nucleus
    make_uniform_soma!(vol_out, 1f0, 0f0)

    # ── Optics ────────────────────────────────────────────────────────────────
    @info "Computing optical propagation…"
    PSF_struct = simulate_optical_propagation(vol_params, psf_params, vol_out)

    # ── Synthetic calcium transient activity ──────────────────────────────────
    @info "Generating synthetic activity ($n_spikes spikes/neuron, prot=$prot)…"
    soma_act, mod_vals, inh_idx = make_naomi_soma_activity(K, nt, dt, n_spikes;
                                                            prot, inh_frac, inh_baseline_boost, seed)
    neur_act   = (soma=soma_act, dend=fill(1f0, K, nt), bg=fill(1f0, K, nt))
    spike_opts = SpikeOpts(; K, nt, dt, prot, N_bg=0, axonflag=false)

    # ── Ground-truth spatial footprints (W_gt) ────────────────────────────────────
    @info "Computing W_gt (individual neuron scans)…"
    sp1  = SpikeOpts(; K, nt=1, dt, prot, N_bg=0, axonflag=false)
    act_base = (soma=fill(1f0,K,1), dend=fill(1f0,K,1), bg=fill(1f0,K,1))
    _, F_base = scan_volume(vol_out, PSF_struct, act_base,
                            ScanParams(; motion=false, verbose=0, scan_buff, sfrac),
                            noise_params, sp1, tpm_params)

    W_gt = zeros(Float64, N1*N2, K)
    for k in 1:K
        soma_k       = fill(1f0, K, 1); soma_k[k,1] = 1f0 + Float32(spike_amp)
        act_k        = (soma=soma_k, dend=fill(1f0,K,1), bg=fill(1f0,K,1))
        _, F_k       = scan_volume(vol_out, PSF_struct, act_k,
                                ScanParams(; motion=false, verbose=0, scan_buff, sfrac),
                                noise_params, sp1, tpm_params)
        W_gt[:, k]   = vec(Float64.(F_k[:,:,1]) .- Float64.(F_base[:,:,1]))
        print("$k ")
    end
    println()

    H_gt = Float64.(soma_act)   # K × nt  (calcium transient traces)

    # ── Normalize W_gt, scale H_gt by spatial power ───────────────────────────────
    # A neuron entirely outside the imaged FOV/focal plane has power 0 — guard
    # the divide so its (already-zero) W_gt column doesn't become NaN.
    powers       = [norm(W_gt[:, k]) for k in 1:K]
    powers_safe  = replace(powers, 0.0 => 1.0)
    W_gt_n  = W_gt ./ powers_safe'
    H_gt_n  = H_gt .* powers

    verbose && @info "W_gt: $(size(W_gt_n))  H_gt: $(size(H_gt_n))"
    verbose && @info "Powers: $(round.(powers; sigdigits=3))"

    vol_out, vol_params, PSF_struct, neur_act, spike_opts, inh_idx, W_gt_n, H_gt_n, powers
end

#== Main scan: generate and save noisy data ==#
function naomi_scan(vol_out, PSF_struct, neur_act, spike_opts, H_gt_n, imgsz, params,
        noise_params, tpm_params; verbose=false)
    @unpack seed, scan_buff, sfrac, pavg, nt = params
    N1, N2 = imgsz

    # NAOMiSim's scan noise (shot/readout) and motion jitter also draw from the
    # global RNG with no seed of their own — reseed here so a noisy scan is
    # reproducible on its own, independent of how many draws prescan consumed.
    seed === nothing ? Random.seed!() : Random.seed!(seed)

    # ── Scan (Poisson shot noise + Gaussian readout noise) ───────────────────────
    scan_params  = ScanParams(; motion=false, verbose=1, scan_buff, sfrac)

    verbose && @info "Scanning volume ($nt frames, pavg=$(pavg) mW…"
    Fnoisy, Fclean = scan_volume(vol_out, PSF_struct, neur_act,
                                scan_params, noise_params, spike_opts, tpm_params)

    _, _, Nt = size(Fclean)
    X_clean = reshape(Float64.(Fclean), N1*N2, Nt) # F is a video tensor (N1, N2, Nt)
    X_noisy = reshape(Float64.(Fnoisy), N1*N2, Nt) # X is a matricized video (unfolded video)

    # SNR: normalize noisy by PMT gain to compare in photon domain
    mu_pmt     = 100.0
    X_noisy_ph = X_noisy ./ mu_pmt
    actual_snr     = 10 * log10(mean(abs2, X_clean) / mean(abs2, X_noisy_ph .- X_clean))
    actual_snr_str = "$(round(actual_snr; digits=1))dB"
    verbose && @info "Actual SNR: $(actual_snr_str)  (pavg=$(pavg) mW)"

    # maxSNR_X: clean-signal frame at peak activity, one column per neuron
    maxindices = [argmax(view(H_gt_n, k, :)) for k in 1:size(H_gt_n, 1)]
    maxSNR_X   = X_clean[:, maxindices]

    X_noisy, X_clean, maxSNR_X, actual_snr, mu_pmt
end