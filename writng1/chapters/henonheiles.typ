= Hénon–Heiles Model

== notes for structre 


- First I define equations of motion and potential and observe trajectories in configuration space.

- notion of perturbed 2D harmonic, with perturbation term $alpha/3 (x^2y+y^3)$ and interaction strength $alpha$ 


$H(p,q)= 1/(2m) (p_x^2+p_y^2)+ (m omega^2)/2 (x^2+y^2)+ alpha/3 (x^2y+y^3)$


- knowing that the phase space captures the dynamics at its fullest, the reduction of the 4D phase space to a 2D surface of section or poincaré section is what I worked on next 

- this should in particular depict the trajectories close to the periodic orbits in configuration space and the surface of section
	- this should make the reader fermiliar with how the information on surface of section relates to the orbits
	- this should lead to the idea of locating periodic orbtis by use of the poincaré map: $M(v)-v=0$
	- here I could also depitct what the contour of $|M(v)-v|$ looks like in proximity of the periodic orbits 

- then we procede with the idea of the monodromy matrix, and to elaborate how we are conisdering function of functions and how 

#pagebreak()



= Hamiltonian 

The Hamiltonian principle is an equivalent formulation to Newtons classical derivation of the equations of motion.
Though is based on powerfull tool, that is used throughout this study – variational principle.
By following the derivation of the Hamiltonian we get insights in this mathematical tool and it may help our understanding.

By Newton the a change in the motion of a particle is the result of a force acting on it. //#cite()

#math.equation(
  alt: "change of motion is induced by a force",
  block: true,
  $m (d^2 bold(x)_i)/(d t^2) = m dot.double(bold(x))_n colon.eq bold(F)_n$,
)


Due to linearity of the forces acting on a body and the force being contiuously differentiable, we can express the force by a potential through a gradient as

#math.equation(
	alt: "force expressed as a gradient",
	block: true,
	$bold(F)_n (bold(x)_n)=- bold(gradient)_n V(bold(x)_n) \
	arrow.l.r.double m dot.double(bold(x))_n+ bold(gradient) V(bold(x)_n)= bold(g)_n (bold(x)_n)=0$
)




The trajectory in configuration space described by $bold(x)(t)$ is the one that sattisfies the Newtons equations of motion.

If we have $N$ particles, then we would have $bold(g)_n (bold(x)_n)$ for $n in [1,dots,N]=I$, and 

#math.equation(
	alt: "",
	block: true,
	$bold(g)((bold(x)_n )_(n in I))=sum_(n in I) bold(g)_n (bold(x)_n)$
)


//#ref("variationalAlgebra")

We then considere the Lagrange–d'Alambert principle of virtual work. We considere variations $bold(delta x)$, so that $bold(delta x)_n (t_1)=bold(delta x)_n (t_2)=0$, for all $n <= N$ particles.

#math.equation(
	alt: "",
	block: true,
	$integral_(t_1) ^t_2 d t space bold(g)(bold(x)) dot bold(delta) bold(x) =&integral_(t_1) ^t_2 d t space m dot.double(bold(x)) dot bold(delta)bold(x)+ integral_(t_1) ^t_2 d t space bold(gradient)V(bold(x)) dot bold(delta)bold(x)\
	=&-integral_(t_1) ^t_2 d t space m (bold(delta) dot(bold(x))^2 )/ 2 + integral_(t_1) ^t_2 d t space bold(delta)V(bold(x)) \
	=& delta integral_(t_1)^(t_2) d t space T-U$
)



#math.equation(
	alt: "",
	block: true,
	$(d dot(bold(x))) / (d t)bold(delta)bold(x)=d / (d t) (dot(bold(x)) bold(delta x))-bold(dot(x) delta) bold(dot(x)) $
)

#math.equation(
	alt: "",
	block: true,
	$bold(delta f(x))=sum_i (partial bold(f(x))) / (partial x_i) bold(delta x_i) $
)

#math.equation(
	alt: "",
	block: true,
	$bold(delta dot(x))^2=2 bold(dot(x) delta x)  $
)

#math.equation(
	alt: "",
	block: true,
	$bold(delta)V(bold(x))=bold(gradient) V(bold(x)) dot bold(delta x)  $
)


So we obtain the well known Lagranian $U((bold(x))_(n in I))-T((bold(dot(x)))_(n in I) ) colon.eq L(bold(dot(x)), bold(x))$ and we can recognize $delta integral_(t_1)^(t_2) d t space L = delta S$ 


- from variational derive $(d) / (d t) (partial L) / (d dot(q))- (partial L)/(partial q)=0$
- with legendre transformation $H(p,q)=p dot(q)-L$


#math.equation(
	alt: "",
	block: true,
	$dot(q) &= (partial H) / (partial p) \
	 dot(p) &= -(partial H) / (partial q)$
)

- this can be reformulated in terms of $bold(u)=(bold(p), bold(q))^T$ such that $dot(bold(u))=bold(f)(bold(u))$

#math.equation(
	alt: "",
	block: true,
	$bold(f)(bold(x))=bold(S)(partial H) / (partial bold(u))$
)

== Periodic Orbits

A periodic orbit is defined by the property, that for a given time $T$ the trajectory $bold(Phi)_(T)(bold(x))=bold(x)$


#figure(
	image("../../figures/henon-heiles/sep16-latex/pot_contour.png", width: 80%),
	caption: [Contour plot of the Hénon–Heiles potential.],
)


#figure(
	image("../../figures/henon-heiles/sep16-latex/orbit_A_eigenvalues.svg", width: 80%),
	caption: [These are the eigenvalues of the monodromy matrix to orbit A, as a funciton of energy.],
)

#figure(
	image("../../figures/henon-heiles/sep16-latex/orbit_A_lyapunov.svg", width: 80%),
	caption: [These are the lyapunov exponets of the monodromy matrix to orbit A.],
)

#figure(
	image("../../figures/henon-heiles/sep16-latex/orbit_A_trace.svg", width: 80%),
	caption: [These are the lyapunov exponets of the monodromy matrix to orbit A.],
)

#figure(
	image("../../figures/henon-heiles/sep16-latex/contour_at_bif_order4.png", width: 80%),
	caption: [Contour plot at the fourth-order bifurcation.],
)

== 