# network_3d: python tools to convert a 3d skeleton to a graph and visualize the result

After finding out that Scikit-image has now a fast [skeletonization algorithm](http://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.skeletonize_3d) I like to share some code to convert such a skeleton to a [networkx](https://networkx.github.io/) graph and a function to visualize the graph in space using mayavi. While the former is strongly influenced by a [matlab code by Philip Kollmannsberger](https://github.com/phi-max/skel2graph3d-matlab) the later is based on a script by [Gael Varoquaux](https://mail.enthought.com/pipermail/enthought-dev/2011-November/030194.html).
Networkx allows does not only allow to quantify the network topology, but also a handy datastructure if you want to write your own analysis tools

## Instalation
I recommend using Anaconda Python and installing current version of 
+ scikit-image
+ networkx
+ mayavi
 


## Usage
At least for the visualization it is best to run `ipython --gui=qt`  
Have a look at the example.py. Without the generation of the data the remaining code usage is straight forward

```python
from skimage.morphology._skeletonize_3d import skeletonize_3d
from network_3d import skel2graph, plot_graph
from mayavi import mlab

#you need to load or create a binary 3D Volume, e.g. see example.py

skel = skeletonize_3d(bin) #scikit-image function
skel = skel.astype(np.bool) #data needs to be bool
G = skel2graph(skel) #create graph

#plot the graph, use the z component to colorcode both the edges and the nodes, scale nodes according to their degree
plot_graph(G,node_color_keyword='z',edge_color_keyword='z',scale_node_keyword='degree')

#show the binary data
mlab.contour3d(bin.astype(np.float),contours=[.5],opacity=.5,color=(1,1,1))
```
This is the result of the visualization:
![Image of 3d graph](https://github.com/refelix/network_3d/blob/master/graph.png?raw=true "result of the example")

 In this example the 'z'-component was used to colorcode the edges and nodes, the size of the nodes is scaled depending on the degree of the nodes.
Have a look at the `plot_graph' docstring for more options: The edge radius can further be varied dependent on a quantity that is stored in the attribute dictionary of the edges.
If you develop your own functions to quantify the edges and you want to visualize this quantity, make sure it is either a float or a list of floats with the same length as the coordinates.

### Interested in more 3D network analysis? 
Currently I am waiting for two scientific publications on the analysis of cellnetworks in bone to be accepted. As soon as this is the case I will make my currently private repository public. In the mean time you can have a look at [my thesis](edoc.hu-berlin.de/docviews/abstract.php?lang=ger&id=42041). If you use this code for scientific work, I would be delighted if you check if my publications are already accepted and acknowledge them in your own work.

