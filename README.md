In contrast to classical configurations (e.g. of structured or unstructured grid) dual resolution meshes consist of a set of superelements with a very simple internal structure. These super-elements are called patches in DrNUM’s internal terminology. Patches can be placed in freely overlapping arrangements and they constitute the actual computational grid. Dual resolution grids can be regarded as unstructured grids of superelements (coarse resolution layer); each superelement, or patch, consists of a highly resolved grid portion which represents the real numerical mesh. At the moment only Cartesian superelements are supported.
 
The goal of this development is to achieve an optimal combination of high numerical and algorithmic efficiency within the superelements, while keeping geometric flexibility on a higher program layer.

This method is thus well suited to run on modern, massively parallel computing devices (e.g. GPUs). It allows to process flow simulations in the order of dozens of millions of cells on classical single node desktop computers. A first example can be seen in the following little video which shows the instability of a compressible jet flow. The grid consists of approximately 16 million cells, fitting into a single 1.5 GB NVidia GTX-580 graphics card; the run time of the case has been approximately 10 hours on one single desktop PC.

http://www.youtube.com/watch?v=S8MBjVtqahY

DrNUM is being developed by numrax GmbH and enGits GmbH.

http://www.numrax.de

http://engits.eu

