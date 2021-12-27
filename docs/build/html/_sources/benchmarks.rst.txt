The Gaussian Kernel is given by :math:`K(i, j) = \sigma^2 \delta_{ij}^2 + \exp(-||x_i - x_j||^2)`. For these benchmarks, we take :math:`\sigma = 10` with :math:`x` being set as a sorted random vector :math:`\in (-1, 1)`. Using the ``plotTree`` function of this library, we can look at the rank structure for this matrix. The following diagram is obtained with :math:`N = 10000`, :math:`M = 500` and tolerance :math:`10^{-12}`

.. image:: images/gaussian_rank_structure.svg
   :width: 400

The green blocks are low-rank blocks. Their intensity of colour shows their degree of "low-rankness". Additionally, the rank has been displayed in each of these blocks. The red blocks are full-rank blocks and would have the rank of :math:`M = 500`

Time Taken vs Tolerance
~~~~~~~~~~~~~~~~~~~~~~~

These benchmarks were performed for size of the matrix :math:`N = 1000000`, with the size of the leaf node set to :math:`M = 100`.

Fast MatVec
^^^^^^^^^^^

+----------------+------------+---------+
|Tolerance       | Assembly(s)|MatVec(s)|
+================+============+=========+
|:math:`10^{-2}` |  12.1156   | 6.97508 |
+----------------+------------+---------+
|:math:`10^{-4}` |  22.649    | 6.79404 |
+----------------+------------+---------+
|:math:`10^{-6}` |  94.8494   | 8.81189 |
+----------------+------------+---------+
|:math:`10^{-8}` |  221.221   | 35.9898 |
+----------------+------------+---------+
|:math:`10^{-10}`|  271.837   | 43.0653 |
+----------------+------------+---------+
..
|:math:`10^{-12}`|     |  |
..
+----------------+------------+---------+
..
|:math:`10^{-14}`|     |  |
..
+----------------+------------+---------+

.. image:: images/gaussian1.png
   :width: 600


Time Taken vs Size of Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For these benchmarks, the leaf size was fixed at :math:`M = 100`, with tolerance set to :math:`10^{-12}`

Fast MatVec
^^^^^^^^^^^

..
+-----------------------+------------+------------+
..
|:math:`N`              | Assembly(s)|MatVec(s)   |
..
+=======================+============+============+
..
|:math:`10^{3}`         |  | |
..
+-----------------------+------------+------------+
..
|:math:`5 \times 10^{3}`|  | |
..
+-----------------------+------------+------------+
..
|:math:`10^{4}`         |   |  |
..
+-----------------------+------------+------------+
..
|:math:`5 \times 10^{4}`|    |   |
..
+-----------------------+------------+------------+
..
|:math:`10^{5}`         |    |    |
..
+-----------------------+------------+------------+
|:math:`5 \times 10^{5}`|     |     |
..
+-----------------------+------------+------------+
..
|:math:`10^{6}`         |     |     |
..
+-----------------------+------------+------------+


.. image:: images/gaussian3.png
   :width: 600
