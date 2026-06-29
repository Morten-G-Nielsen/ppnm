# Multi-Dimensional Non-Linear Root Finder
**Course:** Computational Methods for Physicists  
**Author:** Morten Nielsen  
**Project Task:** Implement Broyden's Method to solve non-linear equations in multi-dimensions.

---

## 1. Implementation
I reused my code from my root finder homework, due to its OOP nature, so it i could easily add two new update rules, one for Good Broyden (GB) and one for Bad Broyden (BB). The reason for implementing both methods was to see how much they differ compared to each other, where the formulas for both methods came from the lecture notes.

Both methods work by doing a rand-1 update to the inverse Jacobian $B$ to reduce the need to do the computationally expensive $\mathcal{O}(N^3)$ QR decomposition. The GB method uses
$$\Delta B=\frac{\Delta x-B\Delta f}{\Delta x^TB\Delta f}\Delta x^TB$$

where GB minimizes the change in the inverse matrix in the step direction. This leads to in preserving curvature leading to the it being able to handle more difficult terrain.


The BB methods uses
$$\Delta B=\frac{\Delta x-B\Delta f}{\Delta f^T\Delta f}\Delta f^T$$

where it minimizes the change in the along the difference $\Delta f_k$. This leads to the update depending heavily on the magnituide of the function output.

In general GB outperforms BB in convergence on complex systems. Due to using this rank-1 update Broyden's methods requires only calculating the Jacobian ones and then updates this leads to a time-complexity of $\mathcal{O}(N^2)$ under ideal circumstances.

To account for situations, where the sugested step size becomes too small the methods are able to recalculate the Jacobian matrix for the system and do the following QR decomposition when $\alpha<1/128$. If the calculated Jacobian does not have linearly independent coloumns the methods will default to identity. This ensures that in the worst case the methods would perform similar to Newton's method.

---

## 2. Testing 

### Test 1: Runtime Scaling & Function Evaluations
To see how the methods compare to Newton's method and to each other a simple function chained function was used that could be of arbitary dimension, where the function had the following form
$$f_i(x)=x_i^2-x_i$$
where $i\in[1,2,3,..., N]$ and the endpoint was set to
$$f_N(x)=x_N^2-1$$
The methods were time with the C++ chrono libary because it was the easiest to implement. It was chosen to not use a log-log plot due to the overlap of the two Broyden methods and also because for low dimensions it lead to a mess.
![Scaling Analysis](plots/plot_scaling.svg)
From the plot that both of Broyden's methods outperforms Newton's method and supriesingly GB performs slightly worse. This makes sense if both methods converge in a similar amount of iterations, so GB will do a small amount of extra calculation.

Both methods were also tested in theire function evaluations, where it can be seen that they vastly outperform Newton's method and performs similar to each other.
![Evaluation Analysis](plots/plot_evals.svg)

### Test 2: Convergence Mapping & Stability Basins
Both methods were then tested on a more complex function as a stress test. The function chosen was a non-linear intersection of a circle and a distorted sinusoid.
$$f_1(x,y)=x^2+y^2-4=0, \quad f_2(x,y)=\sin(x+y^2)-0.5=0$$
The test works by giving an initial guess $(X_0,Y_0)\in[-3,3]$ to map out the iteration count for each method before reaching convergense, with a max iteration count set to 100, so if id did not converge it returned a value of 120 to make it clear on the plots. For this test i chose to turn off the ability to recalculate the Jacobian matrix because i wantet to see how the different update methods handles different start conditions.

From the plot we can see that there are many areas where BB is not able to converge on a solution and also many regions where it uses over half of its max iteration count.
![Stress Test](plots/plot_stability_BB.svg)
For the GB we see that there are still areas where it fails to converge, but there are noticibly fewer and when it is able to converge it does it in very few iterations.
![Stress Test2](plots/plot_stability_GB.svg)
If both methods were allowed to recalclate the Jacobian we would see similar convergence for both methods with the primary difference being in funciton evaluations.

### Test 3: Line Search Trajectory Analysis on Rosenbrock Valley
Due to Newton's method giving the exact curvature makes the line-search method less important. This is not the case for quasi-Newton methods where the line-search method can change drastically on the function evaluations.

To test this GB would find the minimum of the Rosenbrock function using both backtracking and quadratic interpolation line-search. The goal was to be able to see how the quadratic interpolation would follow the curvature of the Rosenbrock function, but the backtracking search perfored in a very similar way so there is not much to see.
![Temp Name](plots/plot_trajectory.svg)

### Test 4: High-Dimensional Physics Application — Static Catenary Equilibrium
A possible application of root finding methods within phyics is in solving equlibrium problems. For this i chose to solve for the shape of a hanging spring-mass chain anchored at $(-5,0)$ and $(5,0)$ under uniform gravity. For $N=100$ moving point masses the state space expands to a $200$-dimensional coupled non-linear system.

For every node $i$ the residual net force must sum to zero at equilibrium
$$\vec{F}_{\text{net}, i} = \vec{F}_{\text{spring}, i} + \vec{F}_{\text{spring}, i+1} + \vec{F}_{\text{gravity}, i} = 0$$
From this requirment the following two equations account for each node.
$$f_{x,i}=k(1-\frac{L}{r_{i-1,i}})(x_i-x_{i-1})-k(1-\frac{L}{r_{i,i+1}})(x_{i+1}-x_i)=0$$
$$f_{y,i}=k(1-\frac{L}{r_{i-1,i}})(y_i-y_{i-1})-k(1-\frac{L}{r_{i,i+1}})(y_{i+1}-y_i)-mg=0$$
Where $k$ is the spring force constant, $L$ is the equlibrium length for each spring, $r_{i,j}$ is the length from node $i$ to node $j$, and $mg$ is the gravitional force with $g=9.81$ and $m=1$.

The solved system can be seen on the figure with a residual force norm within numerical presssicion. In this excat case the max iteration count was set to $300$ to ensure convergences.
![Chain](plots/plot_chain.svg)

---

## 3. Final Remarks & Self-Evaluation

This project would not have been able to cover as much without the use of Google Gemini. It was a great help in the debugging of code, generating the code for gnuplot, and as a way to brainstorm ideas.

I would say i have done the project task and maybe a bit more then nessecary. As i have:
* Implementet both the Good Broyden's method and the Bad Broyden's method.
* Tested their scalabillity agaings Newton's method and eachother.
* Stress tested their convergences on a more complex function.
* Used the Good Broyden in a physical problem where it is nessecary to be able to find the root of a non-linear, multi-dimensional problem.

Of possible improvements would be a better test in line-search methods and how Broyden's methods compare to Newton's methods, but i felt it would go a bit out of scope of the given task. Besides that there are most likely possible optimizations in the code, but again it felt as there are more important things to implement then looking for possible cache misses and fixing that.

A possible critique could be how i have made my interface with the root finding. This is my personal preference to work with, so i think of it as a improvement, but i can see some would prefere maybe a something similar to the methods from scipy. Due to theses reasons i would give myself a score of 9.5/10 depending on how much search methods would have improved the project.
