- 1-parametric-surface-1-field.py plots the flow of one given vector field tangent to a given surface

- 2-parametric-surface-2-fields.py plots one integral curve of two given vector fields tangent to a given surface

- 3-frob-sym.py, given two commuting vector fields, plots a leaf of the integrable distribution they span, 
	computing __symbolocally__ the parametrization given by the composition of flows

5-frob-num.py given two commuting vector fields, plots a leaf of the integrable distribution they span, 
	computing __numerically__ the parametrization given by the composition of flows


figures contains an example plot resulting from each of these programs. 



----

First commit. These programs provide different tastes of the interplay between surfaces and vector fields. They belong to the pre-git era when branches and commits were handled manyally :) hence the abundance of similar files. Maybe I'll clean them in the future. 

Now, 1 is the core of the idea "Plot a vector field of arbitrary components tangent to an arbitrary parametric surface. The core of the idea is to plot a parametric surface with numpy meshgrid, provide the 2 field components, use sympy to find the local basis of tangent vector fields (this is the nice step, making sure the vector field is indeed tangent!), quiver some field values, numerically solve for the flow and plot integral curves. 

File 2 does the same with 2 fields; but these version do the very stupid things of first mapping the vector field from parameter space to the manifold, and then solving the flow ODE. In version 9, the ode is solved in parameter space and the corresponding curve is mapped to the  manifold, which makes much easier to plot curves on spheres or tori, for example.
    
    The other files are on different lines: 3 investigates Frobenious theorem, so it starts with vector fields and finds the parametrization of the tangent distribution. 3 does it simbolically and 5 numerically. 4 is a short interlude about pushforwards, not too relevant (see Sympy for that). 6 is another interlude abot rotations as the acion of Lie groups on manifolds.
    
    7 is nice as it introduces a way to plot a graph surface on a domain which is not simply a product of intervals, so a 2-simplex in R3 can be plotted.
    
    8 puts 1 and 7 togethes, plotting the flow of a vector field on a simplex, with the stupid approach of solving the ode on the simplex.
    
    Recall that 9, as said above, improves 1 by solving the ode in parameter space; so 10 puts 8 and 9 togehter, plotting the flow of a vector field on the simplex, solving the ode in parameter space.
    
    Ok I swear I start using git now and organize things in branches and old commits :)


---

11 considers specifically the replicator equation on the simplex. See replicator folder in game theory for details. 

12 does just plot of parametric surface
