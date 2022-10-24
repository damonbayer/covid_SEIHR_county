using Symbolics
using Latexify
using LinearAlgebra
@variables βₙ βₒ ρ γₙ γₒ νₙ νₒ



T = [
    0 0 βₙ*ρ 0;
    0 0 0 βₒ;
    0 0 0 0;
    0 0 0 0]
-T

Σ = [
    -γₙ 0 0 0;
    0 -γₒ 0 0;
    γₙ 0 -νₙ 0;
    0 γₒ 0 -νₒ
]

K_L = -T * inv(Σ)

latexify(K_L)

E = [1 0; 0 1; 0 0; 0 0]

K = transpose(E) * K_L * E
latexify(K)
K

det(K)

r01 = 1/2 * (tr(K) + sqrt(tr(K)^2 - 4 * det(K)))
r02 = tr(K)/2 *  + sqrt(tr(K)^2/4 - det(K))
r03 = (tr(K) + sqrt((tr(K))^2 - 4 *det(K)))/2

latexify(r03)
latexify(simplify(r03))