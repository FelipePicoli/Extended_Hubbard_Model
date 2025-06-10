"""
    correlation_matrix(psi::MPS,
                       Op1::AbstractString,
                       Op2::AbstractString;
                       kwargs...)

    correlation_matrix(psi::MPS,
                       Op1::Matrix{<:Number},
                       Op2::Matrix{<:Number};
                       kwargs...)

Given an MPS psi and two strings denoting
operators (as recognized by the `op` function),
computes the two-point correlation function matrix
C[i,j] = <psi| Op1i Op2j |psi>
using efficient MPS techniques. Returns the matrix C.

# Optional Keyword Arguments

  - `sites = 1:length(psi)`: compute correlations only
     for sites in the given range
  - `ishermitian = false` : if `false`, force independent calculations of the
     matrix elements above and below the diagonal, while if `true` assume they are complex conjugates.

For a correlation matrix of size NxN and an MPS of typical
bond dimension m, the scaling of this algorithm is N^2*m^3.

# Examples

```julia
N = 30
m = 4

s = siteinds("S=1/2", N)
psi = random_mps(s; linkdims=m)
Czz = correlation_matrix(psi, "Sz", "Sz")
Czz = correlation_matrix(psi, [1/2 0; 0 -1/2], [1/2 0; 0 -1/2]) # same as above

s = siteinds("Electron", N; conserve_qns=true)
psi = random_mps(s, n -> isodd(n) ? "Up" : "Dn"; linkdims=m)
Cuu = correlation_matrix(psi, "Cdagup", "Cup"; sites=2:8)
```
"""





function correlation_matrix(
  psi::MPS, _Op1, _Op2; sites=1:length(psi), site_range=nothing, ishermitian=nothing
)
  if !isnothing(site_range)
    @warn "The `site_range` keyword arg. to `correlation_matrix` is deprecated: use the keyword `sites` instead"
    sites = site_range
  end
  if !(sites isa AbstractRange)
    sites = collect(sites)
  end

  start_site = first(sites)
  end_site = last(sites)

  N = length(psi)
  ElT = scalartype(psi)
  s = siteinds(psi)

  Op1 = _Op1 #make copies into which we can insert "F" string operators, and then restore.
  Op2 = _Op2
  onsiteOp = _op_prod(Op1, Op2)
  fermionic1 = has_fermion_string(Op1, s[start_site])
  fermionic2 = has_fermion_string(Op2, s[end_site])
  if fermionic1 != fermionic2
    error(
      "correlation_matrix: Mixed fermionic and bosonic operators are not supported yet."
    )
  end

  # Decide if we need to calculate a non-hermitian corr. matrix, which is roughly double the work.
  is_cm_hermitian = ishermitian
  if isnothing(is_cm_hermitian)
    # Assume correlation matrix is non-hermitian
    is_cm_hermitian = false
    O1 = op(Op1, s, start_site)
    O2 = op(Op2, s, start_site)
    O1 /= norm(O1)
    O2 /= norm(O2)
    #We need to decide if O1 ∝ O2 or O1 ∝ O2^dagger allowing for some round off errors.
    eps = 1e-10
    is_op_proportional = norm(O1 - O2) < eps
    is_op_hermitian = norm(O1 - dag(swapprime(O2, 0, 1))) < eps
    if is_op_proportional || is_op_hermitian
      is_cm_hermitian = true
    end
    # finally if they are both fermionic and proportional then the corr matrix will
    # be anti symmetric insterad of Hermitian. Handle things like <C_i*C_j>
    # at this point we know fermionic2=fermionic1, but we put them both in the if
    # to clarify the meaning of what we are doing.
    if is_op_proportional && fermionic1 && fermionic2
      is_cm_hermitian = false
    end
  end

  psi = orthogonalize(psi, start_site)
  norm2_psi = norm(psi[start_site])^2

  # Nb = size of block of correlation matrix
  Nb = length(sites)

  C = zeros(ElT, Nb, Nb)

  if start_site == 1
    L = ITensor(1.0)
  else
    lind = commonind(psi[start_site], psi[start_site - 1])
    L = delta(dag(lind), lind')
  end
  pL = start_site - 1

  for (ni, i) in enumerate(sites[1:(end - 1)])
    while pL < i - 1
      pL += 1
      sᵢ = siteind(psi, pL)
      L = (L * psi[pL]) * prime(dag(psi[pL]), !sᵢ)
    end

    Li = L * psi[i]

    # Get j == i diagonal correlations
    rind = commonind(psi[i], psi[i + 1])
    oᵢ = adapt(datatype(Li), op(onsiteOp, s, i))
    C[ni, ni] = ((Li * oᵢ) * prime(dag(psi[i]), !rind))[] / norm2_psi

    # Get j > i correlations
    if !using_auto_fermion() && fermionic2
      Op1 = "$Op1 * F"
    end

    oᵢ = adapt(datatype(Li), op(Op1, s, i))

    Li12 = (dag(psi[i])' * oᵢ) * Li
    pL12 = i

    for (n, j) in enumerate(sites[(ni + 1):end])
      nj = ni + n

      while pL12 < j - 1
        pL12 += 1
        if !using_auto_fermion() && fermionic2
          oᵢ = adapt(datatype(psi[pL12]), op("F", s[pL12]))
          Li12 *= (oᵢ * dag(psi[pL12])')
        else
          sᵢ = siteind(psi, pL12)
          Li12 *= prime(dag(psi[pL12]), !sᵢ)
        end
        Li12 *= psi[pL12]
      end

      lind = commonind(psi[j], Li12)
      Li12 *= psi[j]

      oⱼ = adapt(datatype(Li12), op(Op2, s, j))
      sⱼ = siteind(psi, j)
      val = (Li12 * oⱼ) * prime(dag(psi[j]), (sⱼ, lind))

      # XXX: This gives a different fermion sign with
      # ITensors.enable_auto_fermion()
      # val = prime(dag(psi[j]), (sⱼ, lind)) * (oⱼ * Li12)

      C[ni, nj] = scalar(val) / norm2_psi
      if is_cm_hermitian
        C[nj, ni] = conj(C[ni, nj])
      end

      pL12 += 1
      if !using_auto_fermion() && fermionic2
        oᵢ = adapt(datatype(psi[pL12]), op("F", s[pL12]))
        Li12 *= (oᵢ * dag(psi[pL12])')
      else
        sᵢ = siteind(psi, pL12)
        Li12 *= prime(dag(psi[pL12]), !sᵢ)
      end
      @assert pL12 == j
    end #for j
    Op1 = _Op1 #"Restore Op1 with no Fs"

    if !is_cm_hermitian #If isHermitian=false the we must calculate the below diag elements explicitly.

      #  Get j < i correlations by swapping the operators
      if !using_auto_fermion() && fermionic1
        Op2 = "$Op2 * F"
      end
      oᵢ = adapt(datatype(psi[i]), op(Op2, s, i))
      Li21 = (Li * oᵢ) * dag(psi[i])'
      pL21 = i
      if !using_auto_fermion() && fermionic1
        Li21 = -Li21 #Required because we swapped fermionic ops, instead of sweeping right to left.
      end

      for (n, j) in enumerate(sites[(ni + 1):end])
        nj = ni + n

        while pL21 < j - 1
          pL21 += 1
          if !using_auto_fermion() && fermionic1
            oᵢ = adapt(datatype(psi[pL21]), op("F", s[pL21]))
            Li21 *= oᵢ * dag(psi[pL21])'
          else
            sᵢ = siteind(psi, pL21)
            Li21 *= prime(dag(si[pL21]), !sᵢ)
          end
          Li21 *= psi[pL21]
        end

        lind = commonind(psi[j], Li21)
        Li21 *= psi[j]

        oⱼ = adapt(datatype(psi[j]), op(Op1, s, j))
        sⱼ = siteind(psi, j)
        val = (prime(dag(psi[j]), (sⱼ, lind)) * (oⱼ * Li21))[]
        C[nj, ni] = val / norm2_psi

        pL21 += 1
        if !using_auto_fermion() && fermionic1
          oᵢ = adapt(datatype(psi[pL21]), op("F", s[pL21]))
          Li21 *= (oᵢ * dag(psi[pL21])')
        else
          sᵢ = siteind(psi, pL21)
          Li21 *= prime(dag(psi[pL21]), !sᵢ)
        end
        @assert pL21 == j
      end #for j
      Op2 = _Op2 #"Restore Op2 with no Fs"
    end #if is_cm_hermitian

    pL += 1
    sᵢ = siteind(psi, i)
    L = Li * prime(dag(psi[i]), !sᵢ)
  end #for i

  # Get last diagonal element of C
  i = end_site
  while pL < i - 1
    pL += 1
    sᵢ = siteind(psi, pL)
    L = L * psi[pL] * prime(dag(psi[pL]), !sᵢ)
  end
  lind = commonind(psi[i], psi[i - 1])
  oᵢ = adapt(datatype(psi[i]), op(onsiteOp, s, i))
  sᵢ = siteind(psi, i)
  val = (L * (oᵢ * psi[i]) * prime(dag(psi[i]), (sᵢ, lind)))[]
  C[Nb, Nb] = val / norm2_psi

  return C
end

