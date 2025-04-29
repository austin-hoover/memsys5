# MEMSYS5

This is the source code for the [MEMSYS5 algorithm](https://academic.oup.com/mnras/article/211/1/111/1026907), also known as the "Cambridge Algorithm". The code was downloaded from Steven Gull's [webpage](https://www.mrao.cam.ac.uk/~steve/lisgo/tomography/memsys5/), which no longer exists. I'm unsure if there is another implementation somewhere. I may eventually write a Python interface to the algorithm.

The algorithm uses a conjugate gradient method to search for the maximum-entropy distribution on an $N$-element regular grid. It solves the problem $y = A(x)$ for image vector $x \in \mathbb{R}^N$, given data vector $y \in \mathbb{R}^M$ and forward operator $A$. In most cases $A$ is a sparse matrix. This differs from the [MENT](https://github.com/austin-hoover/ment) algorithm, which uses a continuous formulation of the problem for the specific case of projection data. MENT is the better approach in that case. However, the Cambridge Algorithm is more general and can work for other inverse problems. Another possible advantage is the ability to solve nonlinear inverse problems, such as phase space tomography with strong space charge, and to quantify the reconstruction uncertainty.

**References**:

Skilling, John, and R. K. Bryan. "Maximum entropy image reconstruction: general algorithm." Monthly Notices of the Royal Astronomical Society 211.1 (1984): 111-124.

Myrheim, Jan, and Haavard Rue. "New algorithms for maximum entropy image restoration." CVGIP: Graphical Models and Image Processing 54.3 (1992): 223-238.




