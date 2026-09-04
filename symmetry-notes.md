# Symmetry reduction of periodic orbits in the Hénon–Heiles system

Companion notes to `symmetries.jl`. Each section names the identifiers it
corresponds to in the code.

---

## 1. The system and its point group

The Hamiltonian is

$$
H(x,y,p_x,p_y) \;=\; \frac{p_x^2+p_y^2}{2m} \;+\; V(x,y),
\qquad
V \;=\; \tfrac{1}{2}m\omega^2\!\left(x^2+y^2\right) + a\!\left(x^2y - \tfrac{y^3}{3}\right).
$$

In polar coordinates $x=r\cos\theta$, $y=r\sin\theta$ the cubic term collapses
to a single harmonic. Using $\cos^2\theta\sin\theta-\tfrac13\sin^3\theta
=\sin\theta-\tfrac43\sin^3\theta=\tfrac13\!\left(3\sin\theta-4\sin^3\theta\right)$
and the triple-angle identity $\sin 3\theta = 3\sin\theta-4\sin^3\theta$,

$$
\boxed{\;V(r,\theta) \;=\; \tfrac{1}{2}m\omega^2 r^2 \;+\; \frac{a}{3}\,r^3\sin 3\theta\;}
$$

This single line determines the entire point group.

**Rotations.** $R_\alpha:\theta\mapsto\theta+\alpha$ leaves $V$ invariant iff
$\sin\bigl(3\theta+3\alpha\bigr)=\sin 3\theta$ for all $\theta$, i.e.
$3\alpha\equiv 0 \pmod{2\pi}$, so

$$
\alpha \in \left\{0,\ \tfrac{2\pi}{3},\ \tfrac{4\pi}{3}\right\}
\qquad\Longrightarrow\qquad C_3 .
$$

**Reflections.** A reflection in the line at angle $\varphi$ acts as
$\theta\mapsto 2\varphi-\theta$. Invariance requires
$\sin(6\varphi-3\theta)=\sin 3\theta$ for all $\theta$, hence
$6\varphi-3\theta=\pi-3\theta \pmod{2\pi}$, giving

$$
\varphi \in \left\{\tfrac{\pi}{6},\ \tfrac{\pi}{2},\ \tfrac{5\pi}{6}\right\} .
$$

The choice $\varphi=\pi/2$ is the $y$-axis, i.e. $\sigma: x\mapsto -x$. Note
that $y\mapsto-y$ alone corresponds to $\varphi = 0$, which is **not** in the
list: the cubic term is odd in $y$ at fixed $x$.

Together,

$$
D_3 \;=\; \langle\, C_3,\ \sigma \,\rangle,
\qquad
C_3=\begin{pmatrix}\cos\frac{2\pi}{3} & -\sin\frac{2\pi}{3}\\[2pt] \sin\frac{2\pi}{3} & \cos\frac{2\pi}{3}\end{pmatrix},
\qquad
\sigma=\begin{pmatrix}-1&0\\0&1\end{pmatrix},
\qquad |D_3| = 6 .
$$

> **Code.** `SYM_GROUP` is built from exactly these two generators. The
> statement "$V$ is invariant" is not assumed: `check_group` samples $200$
> random points per element and asserts
> $\bigl|V(Mq)-V(q)\bigr|<10^{-11}$, checks $M^{\!\top}M=\mathbb{1}$, and
> verifies closure by looking up every product $gh$ in the table.

---

## 2. Lifting to phase space, and time reversal

For $M\in O(2)$ define

$$
\Phi_M(q,p) \;=\; (Mq,\;Mp).
$$

This is symplectic: with $\omega=\mathrm{d}p^{\!\top}\!\wedge\mathrm{d}q$,

$$
\Phi_M^{*}\omega \;=\; \mathrm{d}p^{\!\top}M^{\!\top}\wedge M\,\mathrm{d}q
\;=\; \mathrm{d}p^{\!\top}\wedge\mathrm{d}q \;=\; \omega ,
$$

using $M^{\!\top}M=\mathbb{1}$. Since the kinetic term depends only on
$\|p\|$, $H\circ\Phi_M = H$ whenever $V\circ M = V$.

**Time reversal** is

$$
\Theta(q,p) = (q,-p),
$$

which is *anti*-symplectic ($\Theta^*\omega=-\omega$) but still satisfies
$H\circ\Theta=H$; if $t\mapsto u(t)$ solves Hamilton's equations then so does
$t\mapsto \Theta u(-t)$. Because $\Theta$ acts only on the momentum sign it
commutes with every $\Phi_M$, so the full group is a **direct** product:

$$
\boxed{\;G \;=\; D_3\times\mathbb{Z}_2,\qquad |G| = 12\;}
$$

An element is the pair $g=(M,\tau)$ with $\tau=\pm1$, acting as

$$
g\cdot(q,p) \;=\; \bigl(Mq,\ \tau Mp\bigr),
\qquad
(M_1,\tau_1)(M_2,\tau_2) = (M_1M_2,\ \tau_1\tau_2).
$$

> **Code.** `apply_sym(g, u)` implements the action, `compose(g, h)` the
> product by table lookup (which throws if closure fails). Element names are
> $\{\mathrm{E},C_3,C_3^2,\sigma,\sigma C_3,\sigma C_3^2\}$, each optionally
> suffixed by $\Theta$.

---

## 3. What is preserved along a $G$-orbit of periodic orbits

Let $\Gamma$ be a periodic orbit of energy $E$, prime period $T$, crossing the
section $n$ times, with linearised return map $DT\in SL(2,\mathbb{R})$.

For every $g\in G$:

| quantity | reason |
|---|---|
| $E$ | $H\circ g = H$ |
| $n$ | $g$ maps section crossings to section crossings bijectively |
| $T$ | $\Phi_M$ preserves $t$; $\Theta$ reverses it but $\Gamma$ is closed |
| $\operatorname{tr} DT$ | see below |

For the trace, two cases. If $g$ is symplectic, the return map of $g\Gamma$ is
conjugate to that of $\Gamma$ through the (linearised) transition map between
the two sections, and conjugation preserves the trace. If $g$ involves
$\Theta$, the orbit is traversed backwards and $DT\mapsto DT^{-1}$; but in two
dimensions $\det DT = 1$, so

$$
DT^{-1} = \operatorname{adj}(DT)
= \begin{pmatrix} d & -b \\ -c & a\end{pmatrix}
\quad\text{for}\quad
DT=\begin{pmatrix} a & b \\ c & d\end{pmatrix},
\qquad
\operatorname{tr} DT^{-1} = a+d = \operatorname{tr} DT .
$$

So the whole $G$-orbit carries identical $(E,n,T,\operatorname{tr}DT)$, and
hence identical stability class. **This is the justification for analysing one
representative per class**, and simultaneously supplies a cheap necessary
condition for equivalence.

> **Code.** `same_class` tests $n$, $E$, $T$, $\operatorname{tr}DT$ first;
> a mismatch in any of them *proves* the two orbits are inequivalent, so the
> expensive geometric test is skipped.

---

## 4. The surface of section

Define

$$
\Sigma_E \;=\; \bigl\{\,(x,y,p_x,p_y) \;:\; x=0,\ p_x>0,\ H=E \,\bigr\},
$$

coordinatised by $v=(y,p_y)$. Solving $H=E$ at $x=0$,

$$
p_x^2(y,p_y) \;=\; 2m\bigl(E-V(0,y)\bigr)-p_y^2,
\qquad
V(0,y)=\tfrac12 m\omega^2y^2-\tfrac{a}{3}y^3 ,
$$

and the accessible region is $A_E=\{\,p_x^2>0\,\}$, an open set bounded by the
zero-velocity curve. The lift is

$$
\lambda(v) \;=\; \bigl(\,\varepsilon,\ y,\ +\sqrt{p_x^2(v)},\ p_y\,\bigr),
\qquad \varepsilon = 10^{-9},
$$

the small offset placing the state just *past* the section in the direction of
motion, so no spurious crossing is recorded at $t=0$.

The return map $\mathcal{T}:\Sigma_E\to\Sigma_E$ takes $v$ to the next
crossing. A periodic point of order $n$ is a zero of

$$
F_n(v) \;=\; \mathcal{T}^{\,n}(v)-v,
\qquad
DT \;=\; \mathbb{1} + \partial_v F_n \big|_{F_n=0}.
$$

> **Code.** `px2`, `in_section`, `lift`, `poincare_map`, `Fres`, `get_DT`.
> Only upcrossings are recorded (`affect_neg! = nothing` in the
> `ContinuousCallback`), which is exactly the condition $p_x>0$ — the
> orientation convention that the whole coset argument below depends on.

---

## 5. Which symmetries preserve the section?

This is the key structural question, because it decides how much integration
the symmetry reduction costs.

Let $g=(M,\tau)$. Then $g\,\Sigma_E=\Sigma_E$ requires two things.

**(i) $M$ must fix the $y$-axis.** The image of $(0,y)$ must again have zero
$x$-component for all $y$, i.e. $Me_2=\pm e_2$. Since $-\mathbb{1}\notin D_3$
(it is a rotation by $\pi$, not a multiple of $2\pi/3$), only $Me_2=+e_2$ is
possible, which holds for

$$
M \in \{\mathbb{1},\ \sigma\}.
$$

**(ii) $p_x$ must stay positive.** From $p'=\tau Mp$ we get
$p_x' = \tau M_{11}p_x$, and $M_{11}=+1$ for $\mathbb{1}$, $M_{11}=-1$ for
$\sigma$. Positivity therefore forces

$$
(\mathbb{1},+1) = e
\qquad\text{and}\qquad
(\sigma,-1) = \sigma\Theta .
$$

Hence the **section stabiliser** is

$$
\boxed{\;H \;=\; \operatorname{Stab}_G(\Sigma_E) \;=\; \{e,\ \sigma\Theta\},\qquad |H|=2\;}
$$

and its non-trivial element acts on section coordinates as

$$
\sigma\Theta:\ (y,p_y)\;\longmapsto\;(y,-p_y).
$$

Note that this map preserves $p_x^2$, so it maps $A_E$ to itself — the partner
point is automatically admissible, no boundary check needed.

The index is

$$
[G:H] \;=\; \frac{|G|}{|H|} \;=\; \frac{12}{2} \;=\; 6 ,
$$

so $G$ decomposes into **six cosets** $g_iH$, $i=1,\dots,6$. Within a coset the
two elements $g$ and $g\,\sigma\Theta$ give section points related by
$(y,p_y)\mapsto(y,-p_y)$, which costs nothing.

> **Code.** `is_section_stabiliser` implements (i) and (ii) as a predicate and
> `SEC_STAB = filter(is_section_stabiliser, SYM_GROUP)` derives $H$ rather than
> hard-coding it; `COSET_REPS` picks one $g_i$ per coset. The asserts
> `length(SEC_STAB) == 2` and `length(COSET_REPS) == 6` are the numerical
> statement of $|H|=2$, $[G:H]=6$.

---

## 6. Computing the class: five integrations, not eleven

For $g\notin H$ the image $g\,\lambda(v)$ is off the section. Let $\phi^t$ be
the Hamiltonian flow and define the first-return time

$$
t_g \;=\; \inf\bigl\{\,t>t_{\min} \;:\; x\bigl(\phi^t(g\lambda(v))\bigr)=0,\ \dot{x}>0 \,\bigr\},
\qquad
w_g \;=\; \pi\!\left(\phi^{t_g}\!\left(g\lambda(v)\right)\right),
$$

with $\pi(x,y,p_x,p_y)=(y,p_y)$. Because $g\Gamma$ is again a periodic orbit
and $\Sigma_E$ is transversal to the flow, $w_g\in\Sigma_E\cap g\Gamma$: it is
a genuine section point of the image orbit, though not a distinguished one.

The **class root set** is

$$
\mathcal{R}(\Gamma) \;=\; \bigl\{\, w_g \;:\; g\in G \,\bigr\} \subset \Sigma_E,
\qquad |\mathcal{R}| \le 12 .
$$

Cost accounting:

$$
\underbrace{12}_{|G|}
\;\longrightarrow\;
\underbrace{6}_{[G:H]}
\;\longrightarrow\;
\underbrace{5}_{\text{identity coset is free}}
\quad\text{flows per orbit,}
$$

since $w_e=v$ is already known, and each coset's second member follows
algebraically.

> **Code.** `symmetry_data` returns `roots` $=\mathcal{R}(\Gamma)$.
> `flow_to_section` computes $w_g$, discarding crossings with $t<t_{\min}$
> ($10^{-6}$): a symmetry image can be launched arbitrarily close to the
> section, and without this guard the immediate crossing would be picked up
> instead of a meaningful one.

---

## 7. Stabiliser, multiplicity, and orbit–stabiliser

The stabiliser of the *orbit* (not of the section) is

$$
\operatorname{Stab}_G(\Gamma) \;=\; \{\, g\in G \;:\; g\Gamma=\Gamma \,\}.
$$

It is computed by a membership test: since $w_g$ lies on $g\Gamma$, and two
periodic orbits either coincide or are disjoint,

$$
g\Gamma=\Gamma
\quad\Longleftrightarrow\quad
w_g \in \Sigma_E\cap\Gamma
\quad\Longleftrightarrow\quad
\min_{i}\bigl\|w_g - v_i\bigr\| < \mathrm{tol},
$$

where $\{v_i\}=\Sigma_E\cap\Gamma$ are the recorded section points of $\Gamma$.

The **orbit–stabiliser theorem** then gives the number of distinct orbits in
the symmetry class:

$$
\boxed{\;
\mu(\Gamma) \;=\; \bigl|G\cdot\Gamma\bigr| \;=\; \frac{|G|}{\bigl|\operatorname{Stab}_G(\Gamma)\bigr|}
\;=\; \frac{12}{|\operatorname{Stab}_G(\Gamma)|}
\;}
$$

By Lagrange, $|\operatorname{Stab}_G(\Gamma)|$ must **divide** $12$, so

$$
\mu \in \{1,\,2,\,3,\,4,\,6,\,12\}.
$$

This divisibility is a genuine consistency check, not decoration: a stabiliser
of order $5$, say, is impossible and can only mean that the comparison
tolerance is too loose (spurious matches) or too tight (missed matches), or
that a $w_g$ was mis-identified near the energy boundary.

> **Code.** `symmetry_data` accumulates `stab` and computes `mult` via
> `divrem(12, length(stab))`, warning when the remainder is non-zero and
> setting `mult = 0` when some $w_g$ could not be found. `report` then lists
> those rows separately. `sym_type` names the stabiliser subgroup
> (`"C3 invariant"`, `"mirror invariant"`, `"reversible"`, …).

**Worked expectation.** The librating orbit along the $y$-axis lies inside the
mirror plane $x=0$ and retraces itself, so it is fixed by $\sigma$, by
$\Theta$, and by their product — giving $|\mathrm{Stab}|\ge 4$ and $\mu\le3$,
the three rotated copies of the same libration. A generic asymmetric orbit has
trivial stabiliser and $\mu=12$.

---

## 8. The equivalence relation used for deduplication

Define $\Gamma_1\sim\Gamma_2$ iff $\exists g\in G:\ g\Gamma_1=\Gamma_2$. This
is an equivalence relation because $G$ is a group (reflexivity from $e$,
symmetry from $g^{-1}$, transitivity from closure). Operationally,

$$
\Gamma_1\sim\Gamma_2
\quad\Longleftrightarrow\quad
\mathcal{R}(\Gamma_1)\cap\bigl(\Sigma_E\cap\Gamma_2\bigr)\neq\varnothing .
$$

The test is applied only after the invariant prefilter of §3:

$$
n_1=n_2,\quad
E_1=E_2,\quad
T_1\approx T_2,\quad
\operatorname{tr}DT_1\approx\operatorname{tr}DT_2 .
$$

Since $e\in G$ we have $v_1\in\mathcal{R}(\Gamma_1)$, so literal duplicates are
caught by the same test — no separate case needed.

**Restriction to fixed energy.** All comparisons are made within one $E$. Two
section points at different energies may coincide numerically without the
orbits being related, since $\Sigma_{E_1}$ and $\Sigma_{E_2}$ are different
manifolds; only the $\mu$-fold degeneracy at *fixed* $E$ is a symmetry.

The analysis then lives on the quotient

$$
\mathcal{P}_E/G,
$$

where $\mathcal{P}_E$ is the set of periodic orbits at energy $E$. One row per
class is a section of the quotient map.

> **Code.** `same_class` (the relation), `dedup_by` (generic first-wins
> quotient with a configurable predicate), `dedup_classes` (stamps `sym_id`,
> the class label), `dedup_all` (per-energy grouping).

---

## 9. Cost of the reduction, and why the staging matters

Let $N_{\text{seed}}$ be the number of seeds at one energy, $N_{\text{raw}}$
the number of converged rows, $N_{\text{orb}}$ the distinct orbits, and
$N_{\text{cls}}$ the classes. Typically

$$
N_{\text{seed}} \;\gg\; N_{\text{raw}} \;>\; N_{\text{orb}} \;\ge\; \mu\,N_{\text{cls}} .
$$

Annotating symmetry costs $5$ dense integrations per row, so applying it before
deduplication would cost $5N_{\text{raw}}$ flows where $5N_{\text{orb}}$
suffice. The pipeline is therefore staged:

$$
\text{sweep} \;\xrightarrow{\ \text{cheap, } O(N_{\text{raw}}^2)\ \text{point tests}\ }\;
\text{dedup\_raw} \;\xrightarrow{\ 5N_{\text{orb}}\ \text{flows}\ }\;
\text{annotate} \;\xrightarrow{\ \text{invariants first}\ }\;
\text{dedup\_classes}
$$

**Continuation.** Following $N_{\text{cls}}$ representatives instead of
$N_{\text{orb}}$ orbits down a ladder of $|\mathcal{E}|$ energies changes the
cost from

$$
|\mathcal{E}|\cdot N_{\text{orb}}\cdot c_{\text{solve}}
\qquad\text{to}\qquad
|\mathcal{E}|\cdot N_{\text{cls}}\cdot\bigl(c_{\text{solve}}+5c_{\text{flow}}\bigr),
$$

with $N_{\text{orb}}/N_{\text{cls}}=\langle\mu\rangle$ up to $12$. Since
$c_{\text{flow}}\ll c_{\text{solve}}$ (a solve is many Newton steps, each with
a central-difference Jacobian costing four map evaluations), the speedup is
close to $\langle\mu\rangle$.

The full families are recoverable at any point: `expand_symmetry` polishes each
stored $w_g\in\mathcal{R}(\Gamma)$ with one Newton solve — converging in one or
two steps, since $w_g$ is already accurate — deduplicates the results, and
checks the count against $\mu$.

> **Code.** `sweep_energy(...; reduce_to = :classes | :orbits | :raw)`,
> `follow_orbits`, `expand_symmetry`.

---

## 10. Symbol table

| mathematics | code |
|---|---|
| $G=D_3\times\mathbb{Z}_2$, $|G|=12$ | `SYM_GROUP` |
| $g=(M,\tau)$ | `(; name, M, τ)` |
| $g\cdot(q,p)=(Mq,\tau Mp)$ | `apply_sym` |
| $g_1g_2$ | `compose` |
| $H=\operatorname{Stab}_G(\Sigma_E)=\{e,\sigma\Theta\}$ | `SEC_STAB` |
| coset representatives $g_i$, $[G:H]=6$ | `COSET_REPS` |
| $\sigma\Theta:(y,p_y)\mapsto(y,-p_y)$ | `sec_action` |
| $p_x^2(y,p_y)$, $A_E$ | `px2`, `in_section` |
| $\lambda:\Sigma_E\to\mathbb{R}^4$ | `lift` |
| $\mathcal{T}^n$, $F_n$ | `poincare_map`, `Fres` |
| $DT=\mathbb{1}+\partial_vF_n$ | `get_DT` |
| $w_g$ | `flow_to_section` |
| $\mathcal{R}(\Gamma)$ | `sym_y`, `sym_py` |
| $\operatorname{Stab}_G(\Gamma)$ | `sym_stab` |
| $\mu=12/|\operatorname{Stab}|$ | `sym_mult` |
| class label on $\mathcal{P}_E/G$ | `sym_id` |
| $\Gamma_1\sim\Gamma_2$ | `same_class` |
