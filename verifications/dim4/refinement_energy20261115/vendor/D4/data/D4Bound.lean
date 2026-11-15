import Mathlib

/-!
# Sphere packing in four dimensions: density below 0.629521681653

Single-file export of `lean16/Lean16/D4Bound/` (sections LocalEnergy, Spherical, FarPair,
Ledger, Main). Checked with Lean 4 and Mathlib v4.33.1: no `sorry`, axioms
`propext`, `Classical.choice`, `Quot.sound` only.

The main theorem `D4Bound.d4_density_bound` derives the bound from explicit hypotheses
that stand for the inputs certified outside Lean. See `README.md` in this folder and
Section 7 of the paper for exactly what is and is not proved here.
-/

-- ===== LocalEnergy =====
/-!
# The finite local-energy inequality

For a finite family of centres `p : Fin n → V`, a scalar function `f`, a two-point
kernel `K` and an origin constant `A`, the anchored energy of centre `i` is

  E i = A + 2 ∑_{j ≠ i} K (p j - p i) 0 + ∑_{j ≠ i} ∑_{k ≠ i} K (p j - p i) (p k - p i).

If every anchored energy is at least `β`, every pair term
`f d + 2 K d 0 + K d d` is `≤ 0`, and every triangle sum `Γ` is `≤ 0`, then

  ∑_{i,j} f (p j - p i) ≤ n (f 0 + A - β).

This is the finite inequality behind the density bound. Passing from it to the density
of a packing uses Poisson summation and a periodization argument, which are not
formalized here.
-/

namespace D4Bound

open Finset

variable {V : Type*} [AddCommGroup V]

/-- Anchored energy of centre `i`. -/
noncomputable def energy {n : ℕ} (K : V → V → ℝ) (A : ℝ) (p : Fin n → V) (i : Fin n) : ℝ :=
  A + 2 * ∑ j, (if j ≠ i then K (p j - p i) 0 else 0)
    + ∑ j, ∑ k, (if j ≠ i ∧ k ≠ i then K (p j - p i) (p k - p i) else 0)

/-- Triangle sum of the kernel at the three anchors of a triangle. -/
def gamma {n : ℕ} (K : V → V → ℝ) (p : Fin n → V) (i j k : Fin n) : ℝ :=
  K (p j - p i) (p k - p i) + K (p i - p j) (p k - p j) + K (p i - p k) (p j - p k)

/-- Sum over ordered triples of distinct indices of the kernel anchored at the first. -/
noncomputable def tripleSum {n : ℕ} (K : V → V → ℝ) (p : Fin n → V) : ℝ :=
  ∑ i, ∑ j, ∑ k, (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p j - p i) (p k - p i) else 0)

lemma sum_swap12 {n : ℕ} (F : Fin n → Fin n → Fin n → ℝ) :
    ∑ i, ∑ j, ∑ k, F i j k = ∑ i, ∑ j, ∑ k, F j i k := Finset.sum_comm

lemma sum_swap13 {n : ℕ} (F : Fin n → Fin n → Fin n → ℝ) :
    ∑ i, ∑ j, ∑ k, F i j k = ∑ i, ∑ j, ∑ k, F k j i := by
  calc ∑ i, ∑ j, ∑ k, F i j k = ∑ i, ∑ k, ∑ j, F i j k :=
        Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ k, ∑ i, ∑ j, F i j k := Finset.sum_comm
    _ = ∑ k, ∑ j, ∑ i, F i j k := Finset.sum_congr rfl fun k _ => Finset.sum_comm

/-- Three copies of the triple sum are the sum of the triangle sums. -/
lemma three_tripleSum {n : ℕ} (K : V → V → ℝ) (p : Fin n → V) :
    3 * tripleSum K p =
      ∑ i, ∑ j, ∑ k, (if j ≠ i ∧ k ≠ i ∧ k ≠ j then gamma K p i j k else 0) := by
  have h2 : (∑ i, ∑ j, ∑ k,
      (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p i - p j) (p k - p j) else 0)) = tripleSum K p := by
    rw [sum_swap12]
    unfold tripleSum
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ =>
      Finset.sum_congr rfl fun k _ => ?_
    refine if_congr ?_ rfl rfl
    exact ⟨fun h => ⟨Ne.symm h.1, h.2.2, h.2.1⟩, fun h => ⟨Ne.symm h.1, h.2.2, h.2.1⟩⟩
  have h3 : (∑ i, ∑ j, ∑ k,
      (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p i - p k) (p j - p k) else 0)) = tripleSum K p := by
    rw [sum_swap13]
    rw [Finset.sum_congr rfl fun i _ => Finset.sum_comm]
    unfold tripleSum
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ =>
      Finset.sum_congr rfl fun k _ => ?_
    refine if_congr ?_ rfl rfl
    exact ⟨fun h => ⟨h.2.1.symm, h.2.2.symm, h.1⟩, fun h => ⟨h.2.2, h.1.symm, h.2.1.symm⟩⟩
  have hsplit : (∑ i, ∑ j, ∑ k,
      (if j ≠ i ∧ k ≠ i ∧ k ≠ j then gamma K p i j k else 0)) =
      tripleSum K p
      + (∑ i, ∑ j, ∑ k,
          (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p i - p j) (p k - p j) else 0))
      + (∑ i, ∑ j, ∑ k,
          (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p i - p k) (p j - p k) else 0)) := by
    unfold tripleSum
    simp only [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ =>
      Finset.sum_congr rfl fun k _ => ?_
    unfold gamma
    split_ifs <;> ring
  rw [hsplit, h2, h3]
  ring

/-- Separate the diagonal `k = j` from the anchored double sum. -/
lemma split_diag {n : ℕ} (i j : Fin n) (X : Fin n → ℝ) :
    ∑ k, (if j ≠ i ∧ k ≠ i then X k else 0) =
      (if j ≠ i then X j else 0) + ∑ k, (if j ≠ i ∧ k ≠ i ∧ k ≠ j then X k else 0) := by
  have h : ∀ k, (if j ≠ i ∧ k ≠ i then X k else 0) =
      (if j = k then (if j ≠ i then X j else 0) else 0) +
        (if j ≠ i ∧ k ≠ i ∧ k ≠ j then X k else 0) := by
    intro k
    by_cases hk : j = k
    · subst hk; by_cases hji : j = i <;> simp [hji]
    · have hk' : k ≠ j := fun e => hk e.symm
      simp [hk, hk']
  rw [Finset.sum_congr rfl fun k _ => h k, Finset.sum_add_distrib, Finset.sum_ite_eq]
  simp

/-- Grouping: the sum of all anchored energies. -/
lemma sum_energy {n : ℕ} (K : V → V → ℝ) (A : ℝ) (p : Fin n → V) :
    ∑ i, energy K A p i =
      n * A + ∑ i, ∑ j, (if j ≠ i then 2 * K (p j - p i) 0 + K (p j - p i) (p j - p i) else 0)
        + tripleSum K p := by
  have hE : ∀ i, energy K A p i = A + 2 * ∑ j, (if j ≠ i then K (p j - p i) 0 else 0)
      + ∑ j, (if j ≠ i then K (p j - p i) (p j - p i) else 0)
      + ∑ j, ∑ k, (if j ≠ i ∧ k ≠ i ∧ k ≠ j then K (p j - p i) (p k - p i) else 0) := by
    intro i
    unfold energy
    rw [Finset.sum_congr rfl fun j _ => split_diag i j (fun k => K (p j - p i) (p k - p i)),
      Finset.sum_add_distrib]
    ring
  have hP : ∀ i, ∑ j, (if j ≠ i then 2 * K (p j - p i) 0 + K (p j - p i) (p j - p i) else 0) =
      2 * ∑ j, (if j ≠ i then K (p j - p i) 0 else 0)
      + ∑ j, (if j ≠ i then K (p j - p i) (p j - p i) else 0) := by
    intro i
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    split_ifs <;> ring
  rw [Finset.sum_congr rfl fun i _ => hE i, Finset.sum_congr rfl fun i _ => hP i]
  unfold tripleSum
  simp only [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  ring

/-- The finite local-energy inequality. -/
theorem finite_local_energy {n : ℕ} (f : V → ℝ) (K : V → V → ℝ) (A β : ℝ) (p : Fin n → V)
    (hpair : ∀ i j, j ≠ i → f (p j - p i) + 2 * K (p j - p i) 0 + K (p j - p i) (p j - p i) ≤ 0)
    (hgamma : ∀ i j k, j ≠ i → k ≠ i → k ≠ j → gamma K p i j k ≤ 0)
    (hfloor : ∀ i, β ≤ energy K A p i) :
    ∑ i, ∑ j, f (p j - p i) ≤ n * (f 0 + A - β) := by
  -- the triple sum is nonpositive
  have hT : tripleSum K p ≤ 0 := by
    have h := three_tripleSum K p
    have hle : (∑ i, ∑ j, ∑ k,
        (if j ≠ i ∧ k ≠ i ∧ k ≠ j then gamma K p i j k else 0)) ≤ 0 := by
      apply Finset.sum_nonpos; intro i _
      apply Finset.sum_nonpos; intro j _
      apply Finset.sum_nonpos; intro k _
      split_ifs with hc
      · exact hgamma i j k hc.1 hc.2.1 hc.2.2
      · exact le_rfl
    linarith
  -- the energy sum
  have hE : (n : ℝ) * β ≤ ∑ i, energy K A p i := by
    have := Finset.sum_le_sum (fun i (_ : i ∈ (univ : Finset (Fin n))) => hfloor i)
    simpa [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] using this
  rw [sum_energy] at hE
  -- split the diagonal of the scalar sum
  have hf : ∀ i : Fin n, ∑ j, f (p j - p i) =
      f 0 + ∑ j, (if j ≠ i then f (p j - p i) else 0) := by
    intro i
    rw [← Finset.add_sum_erase _ _ (Finset.mem_univ i),
      ← Finset.add_sum_erase _ _ (Finset.mem_univ i)]
    simp only [sub_self, ne_eq, not_true_eq_false, if_false, zero_add]
    congr 1
    refine Finset.sum_congr rfl fun j hj => ?_
    simp [Finset.ne_of_mem_erase hj]
  simp only [hf, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  -- pair terms
  have hp : ∑ i, ∑ j, (if j ≠ i then f (p j - p i) else 0) ≤
      -∑ i, ∑ j, (if j ≠ i then 2 * K (p j - p i) 0 + K (p j - p i) (p j - p i) else 0) := by
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum; intro i _
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum; intro j _
    split_ifs with h
    · linarith [hpair i j h]
    · simp
  linarith

end D4Bound

-- ===== Spherical =====
/-!
# A quartic sum over unit vectors in `ℝ⁴`

For unit vectors `u 0, …, u (n-1)` in `ℝ⁴` with inner products `t i k`, and `a² ≤ 1/2`,

  ∑_{i,k} (t⁴ - a² t²) ≥ (1/8 - a²/4) n².

This is the positive-type bound used by the far-pair lemma. The usual proof writes
`t⁴ - a² t²` in Chebyshev polynomials `U₀, U₂, U₄`, which are positive definite on `S³`.
Here we use an elementary route: with `G = ∑ uᵢ uᵢᵀ` and `T = ∑ uᵢ^{⊗4}`,
`∑ t² = ‖G‖² ≥ n²/4`, and Cauchy–Schwarz against the symmetrized tensor
`P = Sym(G ⊗ δ)` gives `‖T‖² (48 ‖G‖² + 6 n²) ≥ 36 ‖G‖⁴`.
-/

namespace D4Bound

open Finset

/-- The Kronecker delta. -/
def kd (a b : Fin 4) : ℝ := if a = b then 1 else 0

/-- Four-index tensors on `ℝ⁴`. -/
abbrev Idx := Fin 4 × Fin 4 × Fin 4 × Fin 4

/-- The symmetrized trace tensor of a matrix. -/
def symTr (g : Fin 4 → Fin 4 → ℝ) (x : Idx) : ℝ :=
  g x.1 x.2.1 * kd x.2.2.1 x.2.2.2 + g x.1 x.2.2.1 * kd x.2.1 x.2.2.2 +
  g x.1 x.2.2.2 * kd x.2.1 x.2.2.1 + kd x.1 x.2.1 * g x.2.2.1 x.2.2.2 +
  kd x.1 x.2.2.1 * g x.2.1 x.2.2.2 + kd x.1 x.2.2.2 * g x.2.1 x.2.2.1

/-- The fourth tensor power of a vector. -/
def pow4 (v : Fin 4 → ℝ) (x : Idx) : ℝ := v x.1 * v x.2.1 * v x.2.2.1 * v x.2.2.2

lemma sum_idx (F : Idx → ℝ) :
    ∑ x, F x = ∑ a, ∑ b, ∑ c, ∑ d, F (a, b, c, d) := by
  simp only [Fintype.sum_prod_type]

lemma pow4_inner (v w : Fin 4 → ℝ) :
    ∑ x, pow4 v x * pow4 w x = (∑ a, v a * w a) ^ 4 := by
  rw [sum_idx]
  simp only [pow4, Fin.sum_univ_four]
  ring

set_option maxRecDepth 100000 in
lemma pow4_symTr (v : Fin 4 → ℝ) (g : Fin 4 → Fin 4 → ℝ) :
    ∑ x, pow4 v x * symTr g x =
      6 * (∑ a, ∑ b, v a * v b * g a b) * (∑ c, v c ^ 2) := by
  rw [sum_idx]
  simp only [Fin.sum_univ_four]
  simp only [pow4, symTr, kd, Fin.isValue, Fin.reduceEq, if_true, if_false, mul_one, mul_zero,
    zero_mul, one_mul, add_zero, zero_add]
  ring

set_option maxRecDepth 100000 in
lemma symTr_sq (g : Fin 4 → Fin 4 → ℝ) (hg : ∀ a b, g a b = g b a) :
    ∑ x, symTr g x ^ 2 = 48 * ∑ a, ∑ b, g a b ^ 2 + 6 * (∑ a, g a a) ^ 2 := by
  rw [sum_idx]
  simp only [Fin.sum_univ_four]
  simp only [symTr, kd, Fin.isValue, Fin.reduceEq, if_true, if_false, mul_one, mul_zero,
    zero_mul, one_mul, add_zero, zero_add]
  rw [hg 1 0, hg 2 0, hg 3 0, hg 2 1, hg 3 1, hg 3 2]
  ring

lemma trace_sq_bound (g : Fin 4 → Fin 4 → ℝ) :
    (∑ a, g a a) ^ 2 ≤ 4 * ∑ a, ∑ b, g a b ^ 2 := by
  simp only [Fin.sum_univ_four]
  nlinarith [sq_nonneg (g 0 0 - g 1 1), sq_nonneg (g 0 0 - g 2 2), sq_nonneg (g 0 0 - g 3 3),
    sq_nonneg (g 1 1 - g 2 2), sq_nonneg (g 1 1 - g 3 3), sq_nonneg (g 2 2 - g 3 3),
    sq_nonneg (g 0 1), sq_nonneg (g 0 2), sq_nonneg (g 0 3), sq_nonneg (g 1 0),
    sq_nonneg (g 1 2), sq_nonneg (g 1 3), sq_nonneg (g 2 0), sq_nonneg (g 2 1),
    sq_nonneg (g 2 3), sq_nonneg (g 3 0), sq_nonneg (g 3 1), sq_nonneg (g 3 2)]

section family

variable {n : ℕ} (u : Fin n → Fin 4 → ℝ)

/-- Inner products. -/
def ip (i k : Fin n) : ℝ := ∑ a, u i a * u k a

/-- The frame matrix `∑ uᵢ uᵢᵀ`. -/
def frame (a b : Fin 4) : ℝ := ∑ i, u i a * u i b

lemma frame_symm (a b : Fin 4) : frame u a b = frame u b a := by
  unfold frame; exact Finset.sum_congr rfl fun i _ => mul_comm _ _

lemma swap4 {ι κ : Type*} [Fintype ι] [Fintype κ] (F : ι → ι → κ → κ → ℝ) :
    ∑ i, ∑ k, ∑ a, ∑ b, F i k a b = ∑ a, ∑ b, ∑ i, ∑ k, F i k a b := by
  calc ∑ i, ∑ k, ∑ a, ∑ b, F i k a b = ∑ i, ∑ a, ∑ k, ∑ b, F i k a b :=
        Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ a, ∑ i, ∑ k, ∑ b, F i k a b := Finset.sum_comm
    _ = ∑ a, ∑ i, ∑ b, ∑ k, F i k a b :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ a, ∑ b, ∑ i, ∑ k, F i k a b :=
        Finset.sum_congr rfl fun a _ => Finset.sum_comm

lemma ip_sq_expand (i k : Fin n) :
    ip u i k ^ 2 = ∑ a, ∑ b, u i a * u i b * (u k a * u k b) := by
  unfold ip; rw [sq, Finset.sum_mul_sum]
  exact Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => by ring

lemma sum_ip_sq : ∑ i, ∑ k, ip u i k ^ 2 = ∑ a, ∑ b, frame u a b ^ 2 := by
  have h1 := ip_sq_expand u
  have h2 : ∀ a b, frame u a b ^ 2 = ∑ i, ∑ k, u i a * u i b * (u k a * u k b) := by
    intro a b; unfold frame; rw [sq, Finset.sum_mul_sum]
  simp only [h1, h2]
  exact swap4 _

lemma trace_frame (hu : ∀ i, ∑ a, u i a ^ 2 = 1) : ∑ a, frame u a a = n := by
  unfold frame
  rw [Finset.sum_comm]
  simp [← sq, hu]

/-- `∑ t² ≥ n²/4`. -/
lemma sum_ip_sq_ge (hu : ∀ i, ∑ a, u i a ^ 2 = 1) :
    (n : ℝ) ^ 2 ≤ 4 * ∑ i, ∑ k, ip u i k ^ 2 := by
  rw [sum_ip_sq, ← trace_frame u hu]
  exact trace_sq_bound _

/-- The fourth-power tensor of the family. -/
def tens (x : Idx) : ℝ := ∑ i, pow4 (u i) x

lemma tens_sq : ∑ x, tens u x ^ 2 = ∑ i, ∑ k, ip u i k ^ 4 := by
  unfold tens
  simp only [sq, Finset.sum_mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [pow4_inner]; rfl

lemma tens_symTr (hu : ∀ i, ∑ a, u i a ^ 2 = 1) :
    ∑ x, tens u x * symTr (frame u) x = 6 * ∑ i, ∑ k, ip u i k ^ 2 := by
  unfold tens
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  simp only [pow4_symTr, hu, mul_one]
  rw [← Finset.mul_sum]
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [ip_sq_expand]
  unfold frame
  simp only [Finset.mul_sum]
  calc ∑ a, ∑ b, ∑ k, u i a * u i b * (u k a * u k b)
      = ∑ a, ∑ k, ∑ b, u i a * u i b * (u k a * u k b) :=
        Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ = ∑ k, ∑ a, ∑ b, u i a * u i b * (u k a * u k b) := Finset.sum_comm
    _ = ∑ k, ∑ a, ∑ b, u i a * u i b * (u k a * u k b) := rfl

/-- The Cauchy–Schwarz step: `S₄ (48 S₂ + 6 n²) ≥ 36 S₂²`. -/
lemma quartic_cs (hu : ∀ i, ∑ a, u i a ^ 2 = 1) :
    36 * (∑ i, ∑ k, ip u i k ^ 2) ^ 2 ≤
      (∑ i, ∑ k, ip u i k ^ 4) * (48 * ∑ i, ∑ k, ip u i k ^ 2 + 6 * (n : ℝ) ^ 2) := by
  have hcs := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (tens u) (symTr (frame u))
  rw [tens_symTr u hu, tens_sq, symTr_sq _ (frame_symm u), trace_frame u hu,
    ← sum_ip_sq] at hcs
  nlinarith [hcs]

/-- The quartic bound: `∑_{i,k} (t⁴ - a² t²) ≥ (1/8 - a²/4) n²` for `a² ≤ 1/2`. -/
theorem quartic_sum_ge (hu : ∀ i, ∑ a, u i a ^ 2 = 1) (a2 : ℝ) (_ha0 : 0 ≤ a2)
    (ha : a2 ≤ 1 / 2) :
    (1 / 8 - a2 / 4) * (n : ℝ) ^ 2 ≤ ∑ i, ∑ k, (ip u i k ^ 4 - a2 * ip u i k ^ 2) := by
  have hcs := quartic_cs u hu
  have h2 := sum_ip_sq_ge u hu
  set S2 := ∑ i, ∑ k, ip u i k ^ 2
  set S4 := ∑ i, ∑ k, ip u i k ^ 4
  have hsplit : ∑ i, ∑ k, (ip u i k ^ 4 - a2 * ip u i k ^ 2) = S4 - a2 * S2 := by
    simp only [Finset.sum_sub_distrib, ← Finset.mul_sum]
    rfl
  rw [hsplit]
  have hn : (0 : ℝ) ≤ (n : ℝ) ^ 2 := sq_nonneg _
  have hD : 0 < 48 * S2 + 6 * (n : ℝ) ^ 2 + 1 := by nlinarith
  -- (S4 - a2 S2 - c n²) D ≥ 6 (S2 - n²/4) ((6 - 8 a2) S2 + (1/2 - a2) n²) ≥ 0
  have hkey : 0 ≤ (S4 - a2 * S2 - (1 / 8 - a2 / 4) * (n : ℝ) ^ 2) *
      (48 * S2 + 6 * (n : ℝ) ^ 2) := by
    have hf : 0 ≤ (4 * S2 - (n : ℝ) ^ 2) * ((6 - 8 * a2) * S2 + (1 / 2 - a2) * (n : ℝ) ^ 2) := by
      apply mul_nonneg
      · linarith
      · have : 0 ≤ S2 := by nlinarith
        nlinarith
    nlinarith [hf, hcs]
  by_contra hlt0
  have hlt := not_le.mp hlt0
  have hpos : 0 < 48 * S2 + 6 * (n : ℝ) ^ 2 := by
    rcases Nat.eq_zero_or_pos n with h0 | h0
    · subst h0
      exfalso
      have hS : S4 - a2 * S2 = 0 := by simp [S4, S2]
      have h0' : (1 / 8 - a2 / 4) * ((0 : ℕ) : ℝ) ^ 2 = 0 := by simp
      linarith
    · have : (0 : ℝ) < (n : ℝ) ^ 2 := by positivity
      nlinarith
  nlinarith [mul_neg_of_neg_of_pos (sub_neg.mpr hlt) hpos]

/-- Off the diagonal: `∑_{i ≠ k} (t⁴ - a² t²) ≥ (1/8 - a²/4) n² - (1 - a²) n`. -/
theorem quartic_offdiag_ge (hu : ∀ i, ∑ a, u i a ^ 2 = 1) (a2 : ℝ) (ha0 : 0 ≤ a2)
    (ha : a2 ≤ 1 / 2) :
    (1 / 8 - a2 / 4) * (n : ℝ) ^ 2 - (1 - a2) * n ≤
      ∑ i, ∑ k, (if k ≠ i then ip u i k ^ 4 - a2 * ip u i k ^ 2 else 0) := by
  have h := quartic_sum_ge u hu a2 ha0 ha
  have hdiag : ∀ i, ∑ k, (ip u i k ^ 4 - a2 * ip u i k ^ 2) =
      (1 - a2) + ∑ k, (if k ≠ i then ip u i k ^ 4 - a2 * ip u i k ^ 2 else 0) := by
    intro i
    rw [← Finset.add_sum_erase _ _ (Finset.mem_univ i),
      ← Finset.add_sum_erase _ _ (Finset.mem_univ i)]
    have hii : ip u i i = 1 := by unfold ip; simpa [sq] using hu i
    simp only [hii, ne_eq, not_true_eq_false, if_false]
    have : ∑ k ∈ univ.erase i, (if ¬k = i then ip u i k ^ 4 - a2 * ip u i k ^ 2 else 0) =
        ∑ k ∈ univ.erase i, (ip u i k ^ 4 - a2 * ip u i k ^ 2) :=
      Finset.sum_congr rfl fun k hk => by simp [Finset.ne_of_mem_erase hk]
    rw [this]; ring
  simp only [hdiag, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul] at h
  linarith

end family

end D4Bound

-- ===== FarPair =====
/-!
# The far-pair lemma

Let an anchor have near neighbours `yᵢ = rᵢ uᵢ` with `uᵢ` unit vectors in `ℝ⁴`,
`1 ≤ rᵢ ≤ R₅ = 102597/100000`, and `|yᵢ - y_k| ≥ 1`. Put `a = 133/250` and
`q(t) = t²(t² - a²)`. Let `h i k ≥ 0` be the allocation terms, and suppose that a
certified kernel bound gives `h i k ≥ λ q(t_ik)` whenever the pair is far
(`|yᵢ - y_k|² > 49/16`) and `t_ik < -a`. Then

  ∑_{i ≠ k} h i k ≥ λ Q_N,   Q_N = N(13561 N - 179244)/250000.

The geometric facts used are: separated near neighbours have cosine at most `a`, and a
pair with cosine below `-a` is farther apart than `7/4`.
-/

namespace D4Bound

open Finset

/-- The outermost near radius `R₅`. -/
noncomputable def R5 : ℝ := 102597 / 100000

/-- Separated near neighbours have cosine at most `a = 133/250`. -/
lemma cap_le {r s t : ℝ} (hr : 1 ≤ r) (hr' : r ≤ R5) (hs : 1 ≤ s) (hs' : s ≤ R5)
    (hsep : 1 ≤ r ^ 2 + s ^ 2 - 2 * r * s * t) : t ≤ 133 / 250 := by
  unfold R5 at hr' hs'
  have hrs : 0 < r * s := by positivity
  have key : r ^ 2 + s ^ 2 - 1 ≤ 2 * (133 / 250) * (r * s) := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hr) (sub_nonneg.mpr hs),
      mul_nonneg (sub_nonneg.mpr hr) (sub_nonneg.mpr hs'),
      mul_nonneg (sub_nonneg.mpr hr') (sub_nonneg.mpr hs),
      mul_nonneg (sub_nonneg.mpr hr') (sub_nonneg.mpr hs'),
      mul_nonneg (sub_nonneg.mpr hr) (sub_nonneg.mpr hr'),
      mul_nonneg (sub_nonneg.mpr hs) (sub_nonneg.mpr hs')]
  by_contra h
  push Not at h
  nlinarith [mul_lt_mul_of_pos_left h (by positivity : (0 : ℝ) < 2 * (r * s))]

/-- A pair with cosine below `-a` is farther apart than `7/4`. -/
lemma far_of_lt {r s t : ℝ} (hr : 1 ≤ r) (hs : 1 ≤ s) (ht : t < -(133 / 250)) :
    49 / 16 < r ^ 2 + s ^ 2 - 2 * r * s * t := by
  have hrs : 1 ≤ r * s := by nlinarith
  nlinarith [mul_lt_mul_of_pos_left ht (by positivity : (0 : ℝ) < 2 * (r * s))]

/-- `q(t) = t⁴ - a² t² ≤ 0` when `|t| ≤ a`. -/
lemma q_nonpos {t : ℝ} (h1 : -(133 / 250) ≤ t) (h2 : t ≤ 133 / 250) :
    t ^ 4 - (133 / 250) ^ 2 * t ^ 2 ≤ 0 := by
  have : t ^ 2 ≤ (133 / 250) ^ 2 := by nlinarith
  nlinarith [sq_nonneg t]

/-- The far-pair lemma. -/
theorem far_pair_sum {n : ℕ} (u : Fin n → Fin 4 → ℝ) (hu : ∀ i, ∑ a, u i a ^ 2 = 1)
    (r : Fin n → ℝ) (hr1 : ∀ i, 1 ≤ r i) (hr2 : ∀ i, r i ≤ R5)
    (hsep : ∀ i k, k ≠ i → 1 ≤ r i ^ 2 + r k ^ 2 - 2 * r i * r k * ip u i k)
    (lam : ℝ) (hlam : 0 ≤ lam) (h : Fin n → Fin n → ℝ)
    (hnn : ∀ i k, k ≠ i → 0 ≤ h i k)
    (hfar : ∀ i k, k ≠ i → 49 / 16 < r i ^ 2 + r k ^ 2 - 2 * r i * r k * ip u i k →
      ip u i k < -(133 / 250) →
      lam * (ip u i k ^ 4 - (133 / 250) ^ 2 * ip u i k ^ 2) ≤ h i k) :
    lam * ((n : ℝ) * (13561 * n - 179244) / 250000) ≤
      ∑ i, ∑ k, (if k ≠ i then h i k else 0) := by
  -- termwise: h i k ≥ λ q(t_ik)
  have hterm : ∀ i k, k ≠ i →
      lam * (ip u i k ^ 4 - (133 / 250) ^ 2 * ip u i k ^ 2) ≤ h i k := by
    intro i k hki
    have hcap := cap_le (hr1 i) (hr2 i) (hr1 k) (hr2 k) (hsep i k hki)
    by_cases hlt : ip u i k < -(133 / 250)
    · exact hfar i k hki (far_of_lt (hr1 i) (hr1 k) hlt) hlt
    · push Not at hlt
      have := q_nonpos hlt hcap
      nlinarith [hnn i k hki]
  have hsum : lam * ∑ i, ∑ k,
      (if k ≠ i then ip u i k ^ 4 - (133 / 250) ^ 2 * ip u i k ^ 2 else 0) ≤
      ∑ i, ∑ k, (if k ≠ i then h i k else 0) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum; intro i _
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum; intro k _
    split_ifs with hki
    · exact hterm i k hki
    · simp
  have hq := quartic_offdiag_ge u hu ((133 / 250) ^ 2) (by norm_num) (by norm_num)
  have hQ : (1 / 8 - (133 / 250) ^ 2 / 4) * (n : ℝ) ^ 2 - (1 - (133 / 250) ^ 2) * n =
      (n : ℝ) * (13561 * n - 179244) / 250000 := by ring
  rw [hQ] at hq
  nlinarith [mul_le_mul_of_nonneg_left hq hlam]

/-- With `h ≥ 0`, the sum is also at least `λ max(Q_N, 0)`. -/
theorem far_pair_sum_max {n : ℕ} (u : Fin n → Fin 4 → ℝ) (hu : ∀ i, ∑ a, u i a ^ 2 = 1)
    (r : Fin n → ℝ) (hr1 : ∀ i, 1 ≤ r i) (hr2 : ∀ i, r i ≤ R5)
    (hsep : ∀ i k, k ≠ i → 1 ≤ r i ^ 2 + r k ^ 2 - 2 * r i * r k * ip u i k)
    (lam : ℝ) (hlam : 0 ≤ lam) (h : Fin n → Fin n → ℝ)
    (hnn : ∀ i k, k ≠ i → 0 ≤ h i k)
    (hfar : ∀ i k, k ≠ i → 49 / 16 < r i ^ 2 + r k ^ 2 - 2 * r i * r k * ip u i k →
      ip u i k < -(133 / 250) →
      lam * (ip u i k ^ 4 - (133 / 250) ^ 2 * ip u i k ^ 2) ≤ h i k) :
    lam * max ((n : ℝ) * (13561 * n - 179244) / 250000) 0 ≤
      ∑ i, ∑ k, (if k ≠ i then h i k else 0) := by
  have h1 := far_pair_sum u hu r hr1 hr2 hsep lam hlam h hnn hfar
  have h0 : 0 ≤ ∑ i, ∑ k, (if k ≠ i then h i k else 0) := by
    apply Finset.sum_nonneg; intro i _
    apply Finset.sum_nonneg; intro k _
    split_ifs with hki
    · exact hnn i k hki
    · exact le_rfl
  rcases le_total ((n : ℝ) * (13561 * n - 179244) / 250000) 0 with hle | hle
  · rw [max_eq_right hle]; simpa using h0
  · rw [max_eq_left hle]; exact h1

end D4Bound

-- ===== Ledger =====
/-!
# The exact ledger

The 124 count groups `(N, j)` of the inherited census, with lower bounds for their energy
floors, the raw objective and the final arithmetic. `N` is the number of near neighbours and
`j` the first radius `R_j` that contains them. Each floor here is the exact floor of the
certificate (`FINITE_COUNT_ACCEPTANCE.json`, field `profile[].new_floor`) rounded down to a
multiple of `10⁻²⁴`, and `urawUp` is the exact raw objective rounded up. The script
`check_lean_inputs.py` checks these roundings.
-/

namespace D4Bound

/-- `(N, j, lower bound for the energy floor of group (N, j))`. -/
def groups : List (ℕ × ℕ × ℚ) := [
  (0, 1, 14343268241954420762581/50000000000000000000000),
  (1, 1, 271657693628876728298389/1000000000000000000000000),
  (1, 2, 54131850166509922184379/200000000000000000000000),
  (1, 3, 53953045750210461056767/200000000000000000000000),
  (1, 4, 138144490975625700799003/500000000000000000000000),
  (1, 5, 266524898951450114188461/1000000000000000000000000),
  (2, 1, 25674447995401776522533/100000000000000000000000),
  (2, 2, 127428253148199824663811/500000000000000000000000),
  (2, 3, 989068089977385841601/3906250000000000000000),
  (2, 4, 25942206433693990394323/100000000000000000000000),
  (2, 5, 123548217419709893259341/500000000000000000000000),
  (3, 1, 6053139526311162764449/25000000000000000000000),
  (3, 2, 59864247117143378753579/250000000000000000000000),
  (3, 3, 237173828926498810295201/1000000000000000000000000),
  (3, 4, 243091342331757970967777/1000000000000000000000000),
  (3, 5, 228579829740932416787797/1000000000000000000000000),
  (4, 1, 113900427081048974450897/500000000000000000000000),
  (4, 2, 112230277293503096263747/500000000000000000000000),
  (4, 3, 27710284958231424295673/125000000000000000000000),
  (4, 4, 113648336586820293608581/500000000000000000000000),
  (4, 5, 210974940893922989541323/1000000000000000000000000),
  (5, 1, 213770156520907064742349/1000000000000000000000000),
  (5, 2, 1639586421012755206349/7812500000000000000000),
  (5, 3, 103363320245101756102961/500000000000000000000000),
  (5, 4, 106018957050261368618451/500000000000000000000000),
  (5, 5, 24285203192040811165597/125000000000000000000000),
  (6, 1, 100048231465396115424939/500000000000000000000000),
  (6, 2, 39135273522877584243073/200000000000000000000000),
  (6, 3, 19230676863749014836233/100000000000000000000000),
  (6, 4, 197314922350339405572513/1000000000000000000000000),
  (6, 5, 178499740906077900683671/1000000000000000000000000),
  (7, 1, 46690275507315175950349/250000000000000000000000),
  (7, 2, 181888328999206941481089/1000000000000000000000000),
  (7, 3, 1427380170765170299041/8000000000000000000000),
  (7, 4, 18312755516102557676951/100000000000000000000000),
  (7, 5, 6545165769644488326541/40000000000000000000000),
  (8, 1, 173707012187323732744867/1000000000000000000000000),
  (8, 2, 1053142520512654448471/6250000000000000000000),
  (8, 3, 82536877926303456902411/500000000000000000000000),
  (8, 4, 169475669770516235373409/1000000000000000000000000),
  (8, 5, 149669692779364396309853/1000000000000000000000000),
  (9, 1, 160934050642916302225807/1000000000000000000000000),
  (9, 2, 31103929540155243316739/200000000000000000000000),
  (9, 3, 152260329396307012181937/1000000000000000000000000),
  (9, 4, 156359123416746365929727/1000000000000000000000000),
  (9, 5, 34155310939692362417043/250000000000000000000000),
  (10, 1, 37110518658493349197433/250000000000000000000000),
  (10, 2, 142938719493396440511609/1000000000000000000000000),
  (10, 3, 139982099214681567056987/1000000000000000000000000),
  (10, 4, 71888886668825476491989/500000000000000000000000),
  (10, 5, 124483654417262352783997/1000000000000000000000000),
  (11, 1, 68115470699215000491079/500000000000000000000000),
  (11, 2, 65379937948910184042309/500000000000000000000000),
  (11, 3, 128238922545665562975487/1000000000000000000000000),
  (11, 4, 131731476771164981081681/1000000000000000000000000),
  (11, 5, 28314195498194522550711/250000000000000000000000),
  (12, 1, 621502540871105496743/5000000000000000000000),
  (12, 2, 59491487075991491924119/500000000000000000000000),
  (12, 3, 117030656627193984482953/1000000000000000000000000),
  (12, 4, 120220090955223434768349/1000000000000000000000000),
  (12, 5, 10294048372325164647023/100000000000000000000000),
  (13, 1, 4506025287971267057383/40000000000000000000000),
  (13, 2, 21521574298763854469597/200000000000000000000000),
  (13, 3, 106357158697201816124901/1000000000000000000000000),
  (13, 4, 109243473127761298589499/1000000000000000000000000),
  (13, 5, 9353461684661800613167/100000000000000000000000),
  (14, 1, 101281170711546716785599/1000000000000000000000000),
  (14, 2, 154615080258022749007/1600000000000000000000),
  (14, 3, 11745396239456059869/122070312500000000000),
  (14, 4, 12350185065839194636331/125000000000000000000000),
  (14, 5, 2125975965020303843317/25000000000000000000000),
  (15, 1, 22547995237237801236797/250000000000000000000000),
  (15, 2, 86062492392252805737923/1000000000000000000000000),
  (15, 3, 21653473938598911998577/250000000000000000000000),
  (15, 4, 88893970390015194817311/1000000000000000000000000),
  (15, 5, 9681700777971134227347/125000000000000000000000),
  (16, 1, 79382920149430125464857/1000000000000000000000000),
  (16, 2, 75891930424720019719147/1000000000000000000000000),
  (16, 3, 76665800667329105136251/1000000000000000000000000),
  (16, 4, 15785962182337383456577/200000000000000000000000),
  (16, 5, 35389088476711875467737/500000000000000000000000),
  (17, 1, 34426922775459231442061/500000000000000000000000),
  (17, 2, 33061298248300422309281/500000000000000000000000),
  (17, 3, 16751914890204255693169/250000000000000000000000),
  (17, 4, 34493843470726041320737/500000000000000000000000),
  (17, 5, 65012608027711169628291/1000000000000000000000000),
  (18, 1, 58604614391351201750501/1000000000000000000000000),
  (18, 2, 56754347845830264981683/1000000000000000000000000),
  (18, 3, 7218825466456691984101/125000000000000000000000),
  (18, 4, 59446648248565843463771/1000000000000000000000000),
  (18, 5, 29836033293558680043013/500000000000000000000000),
  (19, 1, 48635083908663326609507/1000000000000000000000000),
  (19, 2, 47787041710343265354027/1000000000000000000000000),
  (19, 3, 48894490417773628982163/1000000000000000000000000),
  (19, 4, 5030655207096318429529/100000000000000000000000),
  (19, 5, 51119588135300159978239/1000000000000000000000000),
  (20, 1, 19472555670394911003329/500000000000000000000000),
  (20, 2, 39220535328074830281109/1000000000000000000000000),
  (20, 3, 2527448553569517915391/62500000000000000000000),
  (20, 4, 41567255646579089681547/1000000000000000000000000),
  (20, 5, 42967909436701524425191/1000000000000000000000000),
  (21, 1, 29534553925665672487469/1000000000000000000000000),
  (21, 2, 30454464282739482706231/1000000000000000000000000),
  (21, 3, 8096130071901123352651/250000000000000000000000),
  (21, 4, 33228616213348544168059/1000000000000000000000000),
  (21, 5, 35216887729256437972397/1000000000000000000000000),
  (22, 1, 1275204306326616412341/62500000000000000000000),
  (22, 2, 21427485342817448329679/1000000000000000000000000),
  (22, 3, 23741841532129262292701/1000000000000000000000000),
  (22, 4, 24528229936660260502443/1000000000000000000000000),
  (22, 5, 1690377752640155756967/62500000000000000000000),
  (23, 1, 1443889188175672110267/125000000000000000000000),
  (23, 2, 633981801575736906391/50000000000000000000000),
  (23, 3, 3793740124535520266801/250000000000000000000000),
  (23, 4, 15713114610006354242493/1000000000000000000000000),
  (23, 5, 4692780143331745771753/250000000000000000000000),
  (24, 1, 1891883485914483747101/500000000000000000000000),
  (24, 2, 85169107917086855843/20000000000000000000000),
  (24, 3, 863825531388503827113/125000000000000000000000),
  (24, 4, 7190842437865912643053/1000000000000000000000000),
  (24, 5, 5405612588234897293207/500000000000000000000000),
  (25, 4, 2006898869530340235867/500000000000000000000000),
  (25, 5, 136426225760021236421/31250000000000000000000),
  (26, 5, 3820716433707000462637/1000000000000000000000000)]

/-- Kernel-box multipliers: `λ₁ = 3/10⁶` on `B₁`, `λ₅ = 1/10⁶` on `B₅ ⊇ B_j`. -/
def lam (j : ℕ) : ℚ := if j = 1 then 3 / 1000000 else 1 / 1000000

/-- `Q_N = N²(1/8 - a²/4) - N(1 - a²)` with `a = 133/250`. -/
def QN (N : ℕ) : ℚ := (N : ℚ) * (13561 * N - 179244) / 250000

/-- Raised universal floor (rounded down), attained at group `(24, 1)`. -/
def aLow : ℚ := 1912939165914483747101/500000000000000000000000

/-- The raw objective `(f(0) + A)/8`, rounded up. -/
def urawUp : ℚ := 630000291443570126983569/1000000000000000000000000

lemma QN_formula (N : ℕ) :
    QN N = (N : ℚ) ^ 2 * (1 / 8 - (133 / 250) ^ 2 / 4) - N * (1 - (133 / 250) ^ 2) := by
  unfold QN; ring

set_option maxRecDepth 100000 in
/-- Every group's raised floor is at least `aLow`. -/
theorem groups_ge : ∀ g ∈ groups, aLow ≤ g.2.2 + lam g.2.1 * max (QN g.1) 0 := by
  decide +kernel

/-- The final arithmetic: `urawUp - aLow/8 - 3/8000000 < 0.629521681653`. -/
theorem final_arith : urawUp - aLow / 8 - 3 / 8000000 < 629521681653 / 1000000000000 := by
  unfold urawUp aLow; norm_num

end D4Bound

-- ===== Main =====
/-!
# The four-dimensional density bound, from its certified inputs

This file proves the reduction from the certified inputs to the bound
`Δ₄ < 0.629521681653`. The inputs are explicit hypotheses:

* `hpair`, `hgamma`: the pair and triangle signs of the allocated kernel `K_R`. They are
  the signs of the V7 certificate, checked by exact and interval arithmetic outside Lean;
  the allocation keeps the pair term and makes every triangle sum zero.
* `hraw`: the raw objective `(f(0) + A)/8` of the V7 certificate is at most `urawUp`.
* `hgroup`: each anchor falls in one of the 124 census groups `(N, j)`; its energy without
  the allocation is at least the group floor plus the broad-graph slack `3/10⁶`, and the
  allocation adds at least `λ_j max(Q_N, 0)`. The first part is the inherited count census;
  the second is `far_pair_sum_max` (proved in `FarPair.lean`) applied with the certified
  kernel-box bound `-K ≥ λ_j q` on `B_j`.
* `transfer`: the Cohn–Elkies passage from the finite inequality to density (Poisson
  summation for the Schwartz function `f` with `f̂ ≥ 0`, `f̂(0) = π²/4`, and periodization).

Everything else, including the finite local-energy inequality, the far-pair lemma and the
exact ledger over all 124 groups, is proved here.
-/

namespace D4Bound

open Finset

/-- Squared Euclidean norm on `ℝ⁴`. -/
def sqn (x : Fin 4 → ℝ) : ℝ := ∑ a, x a ^ 2

/-- The finite-patch inequality with the improved constant. -/
theorem d4_patch_bound (f : (Fin 4 → ℝ) → ℝ) (K : (Fin 4 → ℝ) → (Fin 4 → ℝ) → ℝ) (A : ℝ)
    (hpair : ∀ x, 1 ≤ sqn x → f x + 2 * K x 0 + K x x ≤ 0)
    (hgamma : ∀ x y, 1 ≤ sqn x → 1 ≤ sqn y → 1 ≤ sqn (y - x) →
      K x y + K (-x) (y - x) + K (-y) (x - y) ≤ 0)
    (hraw : f 0 + A ≤ 8 * (urawUp : ℝ))
    {n : ℕ} (p : Fin n → (Fin 4 → ℝ)) (hsep : ∀ i j, j ≠ i → 1 ≤ sqn (p j - p i))
    (E H : Fin n → ℝ) (hdec : ∀ i, energy K A p i = E i + H i)
    (hgroup : ∀ i, ∃ g ∈ groups, (g.2.2 : ℝ) + 3 / 1000000 ≤ E i ∧
      ((lam g.2.1 * max (QN g.1) 0 : ℚ) : ℝ) ≤ H i) :
    ∑ i, ∑ j, f (p j - p i) ≤ n * (8 * (urawUp : ℝ) - aLow - 3 / 1000000) := by
  have hfloor : ∀ i, (aLow : ℝ) + 3 / 1000000 ≤ energy K A p i := by
    intro i
    obtain ⟨g, hg, hE, hH⟩ := hgroup i
    have hq : ((aLow : ℚ) : ℝ) ≤ ((g.2.2 + lam g.2.1 * max (QN g.1) 0 : ℚ) : ℝ) :=
      Rat.cast_le.mpr (groups_ge g hg)
    rw [Rat.cast_add] at hq
    rw [hdec i]
    linarith
  have hg' : ∀ i j k, j ≠ i → k ≠ i → k ≠ j → gamma K p i j k ≤ 0 := by
    intro i j k hji hki hkj
    have h := hgamma (p j - p i) (p k - p i) (hsep i j hji) (hsep i k hki)
      (by have := hsep j k hkj; rwa [show p k - p i - (p j - p i) = p k - p j by abel])
    unfold gamma
    have e1 : -(p j - p i) = p i - p j := by abel
    have e2 : p k - p i - (p j - p i) = p k - p j := by abel
    have e3 : -(p k - p i) = p i - p k := by abel
    have e4 : p j - p i - (p k - p i) = p j - p k := by abel
    rw [e1, e2, e3, e4] at h
    exact h
  have hmain := finite_local_energy f K A ((aLow : ℝ) + 3 / 1000000) p
    (fun i j hji => hpair _ (hsep i j hji)) hg' hfloor
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have : (n : ℝ) * (f 0 + A - ((aLow : ℝ) + 3 / 1000000)) ≤
      n * (8 * (urawUp : ℝ) - aLow - 3 / 1000000) := by
    apply mul_le_mul_of_nonneg_left _ hn; linarith
  linarith

/-- **The bound.** Under the certified inputs, the density of every sphere packing in `ℝ⁴`
is below `0.629521681653`. -/
theorem d4_density_bound (f : (Fin 4 → ℝ) → ℝ) (K : (Fin 4 → ℝ) → (Fin 4 → ℝ) → ℝ) (A : ℝ)
    (density : ℝ)
    (transfer : ∀ C : ℝ, (∀ (n : ℕ) (p : Fin n → (Fin 4 → ℝ)),
      (∀ i j, j ≠ i → 1 ≤ sqn (p j - p i)) → ∑ i, ∑ j, f (p j - p i) ≤ n * C) →
      density ≤ C / 8)
    (hpair : ∀ x, 1 ≤ sqn x → f x + 2 * K x 0 + K x x ≤ 0)
    (hgamma : ∀ x y, 1 ≤ sqn x → 1 ≤ sqn y → 1 ≤ sqn (y - x) →
      K x y + K (-x) (y - x) + K (-y) (x - y) ≤ 0)
    (hraw : f 0 + A ≤ 8 * (urawUp : ℝ))
    (hcensus : ∀ (n : ℕ) (p : Fin n → (Fin 4 → ℝ)), (∀ i j, j ≠ i → 1 ≤ sqn (p j - p i)) →
      ∃ E H : Fin n → ℝ, (∀ i, energy K A p i = E i + H i) ∧
        ∀ i, ∃ g ∈ groups, (g.2.2 : ℝ) + 3 / 1000000 ≤ E i ∧
          ((lam g.2.1 * max (QN g.1) 0 : ℚ) : ℝ) ≤ H i) :
    density < 629521681653 / 1000000000000 := by
  have hd := transfer (8 * (urawUp : ℝ) - aLow - 3 / 1000000) (by
    intro n p hsep
    obtain ⟨E, H, hdec, hgroup⟩ := hcensus n p hsep
    exact d4_patch_bound f K A hpair hgamma hraw p hsep E H hdec hgroup)
  have hfin : ((urawUp - aLow / 8 - 3 / 8000000 : ℚ) : ℝ) <
      ((629521681653 / 1000000000000 : ℚ) : ℝ) := Rat.cast_lt.mpr final_arith
  push_cast at hfin
  linarith

end D4Bound
