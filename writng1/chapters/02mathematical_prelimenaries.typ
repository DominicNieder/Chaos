= Mathematical prelimenaries

- introduce vector spaces
- differential calculus
  - $delta f(x)= d/(d epsilon) f(x+epsilon eta)$ #label("variationalAlgebra")
  - fundamental lemma of variational calculus
- dynamical systems
- topology

- symplectic algebra
...

= dynamical systems

A dynamical system $(Phi_t)_(t in G)$ is defined over a phasespace $M$ with elements $bold(x) in M$. Chose a time $t$ and a some starting point and the dynamical system will defines where you end up:



#math.equation(
	alt: "continuous dynamical system",
	block: true,
	$bold(Phi): G times M arrow.r.long M, space.quad bold(Phi) (t,bold(x)) colon.eq bold(Phi)_t (bold(x))$
)

To be precise, $G$ is a topological group, like $( RR, +)$ and to permit a continuouse dynamical system, $M$ is required to be a Hausdorff-space, and $G times M$ makes a topological product. 

In the study of dynamical systems we are interested in the behaviour of $Phi$ that is governed by

#math.equation(
	alt: "continuous dynamical system",
	block: true,
	$(d bold(Phi)_t )/(d t) equiv dot(bold(Phi))=bold(f)(bold(x)).$
)
