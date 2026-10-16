= Hénon–Heiles Model

== notes for structre 


- First I define equations of motion and potential and observe trajectories in configuration space.

- notion of perturbed 2D harmonic, with perturbation term $alpha/3 (x^2y+y^3)$ and interaction strength $alpha$ 


$H(p,q)= 1/(2m) (p_x^2+p_y^2)+ (m omega^2)/2 (x^2+y^2)+ alpha/3 (x^2y+y^3)$

#figure(
	image("../../figures/henon-heiles/sep16-latex/pot_contour.png", width: 80%),
	caption: [Contour plot of the Hénon–Heiles potential.],
)

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
  $m (d^2 x)/(d t^2) colon.eq bold(F)$,
)


Due to linearity of the forces acting on a body, we can express the force by a potential through a gradient as

#math.equation(
	alt: "force expressed as a gradient",
	block: true,
	$bold(F)=- bold(gradient) V(bold(x)).$
)

Then by following derivation is based on Helmholz. We can considere the reformulation

#math.equation(
	alt: "",
	block: true,
	$bold(f)(bold(x))=m d^2 / (d t^2) bold(x)+ bold(gradient) V(bold(x))=0$
)

The trajectory in configuration space described by $bold(x)(t)$ is the one that sattisfies the Newtons equations of motion.


In order to find a solution of this one makes use of 
when we consider the Functional $integral_(t_1) ^t_2 bold(f)(bold(x)) d t$
//#ref("variationalAlgebra")

We then considere the Helmholz, the work 


#math.equation(
	alt: "",
	block: true,
	$integral_(t_1) ^t_2 bold(f)(bold(x)) d t=integral_(t_1) ^t_2 $
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