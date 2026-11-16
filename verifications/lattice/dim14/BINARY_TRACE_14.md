# Direct binary obstruction, independent of spectral distortion

For the 14D frame normalization G_ii=2+a_i, |G_ij|<=(a_i+a_j)/2, let S=sum a_i. Assume S<4. Physical squared norm of coordinate vector x is x^T G x/4.

The integral odd-sum coset, if present, contains every e_i modulo D14. Minimum1 would imply G_ii/4>=1 for each i, so a_i>=2 and S>=28. Thus the coordinate-residue map has no kernel.

An order-two element has half-integer support A of even cardinality w. Its representatives on that support have coordinates epsilon_i/2, with signs of prescribed parity, and zero elsewhere. For w=4 or6, averaging over all sign patterns of the required parity annihilates every off-diagonal term. Minimum1 therefore gives

    1 <= (1/16) sum_{i in A} G_ii,
    sum_{i in A} a_i >= 16-2w.

The complement has at least two coordinates, and its total excess is nonnegative: summing its pair inequalities gives (|Ac|-1)sum_{Ac}a_i>=0. Consequently S>=8 for w4 and S>=4 for w6.

For w2, take the applicable representative (e_i +/- e_j)/2. Its physical norm is (G_ii+G_jj +/-2G_ij)/16 <= (4+2(a_i+a_j))/16. Minimum1 implies a_i+a_j>=6, hence S>=6.

Thus when S<4, every nonzero binary word has weight>=8. A rank-four subcode would have fifteen nonzero words of total weight>=120; each of fourteen coordinates contributes at most eight ones, so total<=112. Therefore binary rank<=3. This improves the range of the binary step alone; cyclic order-eight exclusions still require separate metric certificates, and no stronger full14D trace theorem follows until those are checked.
