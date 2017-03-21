
import networkx as nx
import collections
import itertools as IT
from scipy.ndimage import label
from mayavi import mlab
import numpy as np
import sys
import pdb

def get_nhood(img):
    """
    calculates the neighborhood of all voxel, needed to create graph out of skel image
    inspired by [Skel2Graph](https://github.com/phi-max/skel2graph3d-matlab) from Philip Kollmannsberger.
    Parameters
    ----------
    img : ndarray
        Binary input image.  Needs to be padded with zeros
    
    Returns
    -------
    nhood : neighborhood of all voxels that are True
    nhi : the indices of all 27 neighbors, incase they are true
    inds : of all True voxels (raveled)
    """
    print 'calculating neighborhood'
    #check if padded
    assert img.dtype==np.bool
    assert img[0,:,:].max()==0
    assert img[-1,:,:].max()==0
    assert img[:,0,:].max()==0
    assert img[:,-1,:].max()==0
    assert img[:,:,0].max()==0
    assert img[:,:,-1].max()==0

    inds=np.arange(img.size)
    inds=inds[img.ravel()]
    
    x, y, z = np.unravel_index(inds, img.shape)  # unraveled indices
    # space for all neighbors
    nhi = np.zeros((inds.shape[0], 27), dtype=np.uint32)
    
    s0, s1, s2 = np.array(img.shape) -2
    assert s0 * s1 * s2 < np.iinfo(nhi.dtype).max, 'the array is too big, fix datatype of nhi'
        
    
    for xx in range(0, 3):
        for yy in range(0, 3):
            for zz in range(0, 3):
                n = img[xx:(s0 + xx), yy:(s1 + yy), zz:(s2 + zz)]  # shifted image
                
                w = np.lib.ravel_multi_index(np.array([xx, yy, zz]), (3, 3, 3))
                
                nhi[:, w] = (np.lib.ravel_multi_index(np.array([xx + x - 1, yy + y - 1, zz + z - 1]), np.array(img.shape)) * n[img[1:-1, 1:-1, 1:-1]])

                
    # calculate number of neighbors (+1)
    nhood = np.sum(nhi > 0, axis=-1, dtype=np.uint8)
    
    print 'done'
    return nhood, nhi, inds


def skel2graph(skel,  origin=[0, 0, 0], spacing=[1, 1, 1], pad=True):
    """
    converts skeletonized image into networkx structure with coordinates at node and edge attributes.
    This algorithmn is based on the matlab algorithmn of [Philip Kollmannsberger](https://github.com/phi-max/skel2graph3d-matlab) and was written to analyze the denspacingitic cell networks.
nodevoxels (neighborhood >3) which are not seperated by edge voxels (neighborhood = 2) will become one node with coordinates 
    Parameters
    ----------
    skel : ndarray (bool)
        Binary skeleton image.  
    origin : list, tuple
        coordinates of the corner of the 3D image
    spacing : list, tuple
        spacing of the voxels. important to give length of the edges a meanigful interpretation
    pad : bool
        should the skeleton be padded first. If you are unsure if the skeleton was padded before, you should enable this option
    
    
    Returns
    -------
    G : networkx.Graph
        
        """
    
    assert skel.dtype==np.bool
    if pad:
        img = np.pad(skel, pad_width=1, mode='constant')
        
        origin -= np.array(spacing)
    else:
        img = skel
    
    G = nx.MultiGraph()
    w, l, h = np.shape(img)

    # get neighborhood
    nhood, nhi, inds = get_nhood(img)
 
    # edge voxesl are those with 2 neighbors
    e_v = nhood == 3
    e_inds = inds[e_v] #position

    # all others are node voxels
    n_v = (nhood != 3) 
    n_inds = inds[n_v] #position
    
    # what are the two neighbors
    e_nh_idx = nhi[e_v]
    e_nh_idx[:, 13] = 0
    e_nb = np.sort(e_nh_idx)[:, -2:]#find the two neighbors
    #free memory    
    del nhi    
    del e_nh_idx

    # get clusters of node voxels
    nodes = np.zeros_like(img, dtype=bool)
    nodes[np.lib.unravel_index(inds[n_v], [w, l, h])] = 1
    indnode, num = label(nodes, structure=np.ones((3, 3, 3)))
    #make collection of all nodevoxel in cluster
    node_inds_all = collections.defaultdict(list)
    if num < 2**16:
        indnode = indnode.astype(np.uint16)
    try: #this might produce memory error
        inspacingavel = indnode.ravel()

        for idx, val in IT.izip(np.arange(nodes.size, dtype=np.uint64), inspacingavel):
            stop
            if not val == 0:
                node_inds_all[val].append(idx)

        del indnode
    except:  # do it in juncks
        
        inspacingavel = indnode.ravel()
        del indnode
        steps = np.linspace(
            0, inspacingavel.size, inspacingavel.size / 10000, dtype=np.uint64)
        for i in range(len(steps) - 1):
            for idx, val in IT.izip(np.arange(steps[i], steps[i + 1], dtype=np.uint64), inspacingavel[steps[i]:steps[i + 1]]):
                if not val == 0:
                    node_inds_all[val].append(idx)
    N=len(node_inds_all.keys())
    for i in node_inds_all.keys():
        sys.stdout.write("\rpercentage: " + str(100 * i / N) + "%")
        sys.stdout.flush()
        node_inds = node_inds_all[i]
        i -= 1  # to have n_v starting at 0
        x, y, z = np.lib.unravel_index(node_inds, nodes.shape)

        xmean = x.mean()
        ymean = y.mean()
        zmean = z.mean()
       
        if len(x) > 1:
            # node concids out of more than one voxel (only important for
            # debugging?)
            mn = True
        else:
            mn = False
        G.add_node(i, {'x': xmean * spacing[0] + origin[0], 'y': ymean * spacing[1] + origin[1], 'z': zmean * spacing[
                   2] + origin[2], 'multinode': mn,  'idx': node_inds})
        # find all canal vox in nb of all node idx
    print 'nodes done, track edges'
    # del node_inds_all
    nodex = nx.get_node_attributes(G, 'x')
    nodey = nx.get_node_attributes(G, 'y')
    nodez = nx.get_node_attributes(G, 'z')
    edge_inds = []

    # now, find all edges

    e_ravel = e_nb.ravel()  # might make it faster
    
    N=len(G.nodes())
    alledges=[]
    for i, nd in G.nodes(data=True):
        sys.stdout.write("\rpercentage: " + str(100 * i / N) + "%")
        sys.stdout.flush()
        # if i ==46:
        # pdb.set_trace()
        node_inds = nd['idx']
        if nhood[inds==node_inds[0]]==1:
            continue
        # find all edge voxels which are neighbors of node i
        nbs = np.in1d(e_ravel, node_inds).reshape(e_nb.shape)
        e_n, pos = np.where(nbs)
        # iterate through edge untill node is reached
        for m, idx in enumerate(e_n):
            edge_inds = [e_inds[idx]]
            test = e_nb[idx, pos[m] == 0]
            # pdb.set_trace()
            while test in e_inds:
                edge_inds.append(test)
                newcan = e_nb[e_inds == test]
                otest = test

                test = newcan[newcan != edge_inds[-2]][0]
            # pdb.set_trace()
            if test in n_inds:
               
                n2 = inspacingavel[test] - 1
                
                # to avoid starting the edge from the oposite side
                e_nb[e_inds == edge_inds[-1]] *= 0
                assert  n2  in G.nodes()
                
                x, y, z = np.lib.unravel_index(
                    np.array(edge_inds).astype(np.int64), img.shape)
                x_ = np.r_[nodex[i], x * spacing[0] + origin[0], nodex[n2]]
                y_ = np.r_[nodey[i], y * spacing[1] + origin[1], nodey[n2]]
                z_ = np.r_[nodez[i], z * spacing[2] + origin[2], nodez[n2]]
                alledges.append((i, n2, {'x': x_, 'y': y_, 'z': z_, 'weight': get_length(
                    x_, y_, z_), 'n1': i, 'n2': n2}))
              
    G.add_edges_from(alledges)
    print 'edges done'
    

    return G



def get_r(inds, shape, origin=[0, 0, 0], spacing=[1, 1, 1], order='C'):
    """
    get position of voxels 
    Parameters
    ----------
    inds : ndarray
        Binary skeleton image.  
    origin : list, tuple
        coordinates of the corner of the 3D image
    spacing : list, tuple
        spacing of the voxels. important to give length of the edges a meanigful interpretation
    order : bool
        see np.lib.unravel_index documentation
    
    
    Returns
    -------
        r : list of ndarrays
            position of voxels

        """
    return map(np.add, origin, map(np.multiply, spacing, np.lib.unravel_index(inds, shape, order=order)))


def get_cms(inds, shape, r0=[0, 0, 0], dr=[1, 1, 1], order='C'):
    """
    get mean position of voxels
      
    Parameters
    ----------
    inds : ndarray
        Binary skeleton image.  
    origin : list, tuple
        coordinates of the corner of the 3D image
    spacing : list, tuple
        spacing of the voxels. important to give length of the edges a meanigful interpretation
    order : bool
        see np.lib.unravel_index documentation
    
    
    Returns
    -------
        r : tuple 
            mean position of voxels

        """
    x_, y_, z_ = get_r(inds, shape, r0, dr, order=order)
    return x_.mean(), y_.mean(), z_.mean()


def get_length(x_, y_, z_):
    """
    calculate length of an edge
    Parameters:
    -----------
    x_ : list
        x coordinate 
    y_ : list
        y coordinate 
    z_ : list
        z coordinate

    Returns
    -------
    length : float
        length calculated as sum od eucledian distances between voxels in list
        

    """
    if len(x_) > 1:
        return np.sum(np.sqrt((np.array(x_[1:]) - np.array(x_[:-1]))**2 + (np.array(y_[1:]) - np.array(y_[:-1]))**2 + (np.array(z_[1:]) - np.array(z_[:-1]))**2))
    else:
        return 0



def plot_graph( G,
               node_size=1, node_color=(1, 1, 1), scale_node_keyword=None,
               node_color_keyword=None, tube_radius=.4, edge_color=(1, 1, 1),
               edge_color_keyword=None, edge_radius_keyword=None,
                edge_colormap='jet', **kwargs):
    """ 3D plot of a 3D  network, inspired by a code snipped from from Gael Varoquaux posted at
https://mail.enthought.com/pipermail/enthought-dev/2011-November/030194.html
  opposed to the original code, where nodes are connected by straight lines, this function uses a list of coordinates to visualize a network which might represent a 3D skeleton
 For both, edges and nodes the size and the color can be used to visualize a parameter of the attribute dictionary, for edges this needs to be either a a number per edge, or a sequence with the length equal to the length of coordinates

        Parameters
        ----------
        G : networkx.Graph
            nodes and edges must have coordinates stored in attribute dictionary named 'x','y','z'
        node_size : float 
            size of spheres 
        node_color : tuple
            color of sphears
        edge_color : tuple
            color of tubes

        scale_node_keyword : string or None 
            if None, constant sizes are used, otherwise 
                the nodes are scaled with float given in G.node[i][scale_node_keyword] could also be 'degree'
        node_color_keyword: string or None
            if None is given node spheres have the same color, otherwise the float value of G.node[i][node_color_keyword] is used in cobination with a defauld colormap
        edge_color_keyword : string or None
            if None use edgecolor, otherwise G[i][j][edge_color_keyword] is used to colorecode this value or list of values
        edge_radius_keyword : string or None
            if None use edge_radius, otherwise Ge[i][j][edge_radius_keyword] is used to vary the radius according to  this value or list of values
        Returns
        -------
        tube_surf : tvtk actor
            actor of the edges
        pts : tvtk actor
            actor of the nodes
    """ 
    print 'plot graph'
    from tvtk.api import tvtk
    cell_array = tvtk.CellArray()
   
    # initialize coordinates of edge segments
    xstart = []
    ystart = []
    zstart = []
    xstop = []
    ystop = []
    zstop = []
    edge_r = [] #for radius
    edge_c = [] #for color
    edgecount = 0
       
    lines = []
    # initialize coordinates of nodes
    xn = []
    yn = []
    zn = []
  

    for node in G.nodes():
        xn.append(G.node[node]['x'])
        yn.append(G.node[node]['y'])
        zn.append(G.node[node]['z'])
       
    
    # add edge segments
    i = 0 #number of added segments
    for n1, n2, edgedict in G.edges(data=True):

        edgecount += 1
        xstart.extend(edgedict['x'])
        xstop.extend(edgedict['x'])
        ystart.extend(edgedict['y'])
        ystop.extend(edgedict['y'])
        zstart.extend(edgedict['z'])
        zstop.extend(edgedict['z'])
        #how many segments in line?
        l = len(edgedict['x'])
        line = range(i, i + l)
        line.insert(0, l)
        lines.extend(line)
        i += l
        #color
        if edge_color_keyword is None:
            edge_c += [1] * (l)

        elif  edge_color_keyword not in edgedict.keys():
            edge_c += [1] * (l)

        else:
            
            if len(edgedict[edge_color_keyword]) == 1:
                edge_c.extend([edgedict[edge_color_keyword]] * (l))
            
            else:
                assert len(edgedict[edge_color_keyword]) == len(edgedict['x'])
                edge_c.extend(edgedict[edge_color_keyword].tolist())

        #radius
        if edge_radius_keyword in edgedict.keys():
            if np.size(edgedict[edge_radius_keyword]) == 1:
                edge_r.extend([edgedict[edge_radius_keyword]] * (l))
            else:
                assert len(edgedict[edge_radius_keyword]) == len(edgedict['x'])
                edge_r.extend(edgedict[edge_radius_keyword])
        else:
            edge_r += [tube_radius] * (l)
    #pdb.set_trace()
    #pot network
    fig = mlab.gcf()
    disable_render = fig.scene.disable_render
    fig.scene.disable_render = True

    edge_c = np.array(edge_c)
    if np.isinf(edge_c).any():
        edge_c[np.isinf(edge_c)] = np.nan

    

    xv = np.array([xstart])  
    yv = np.array([ystart])  
    zv = np.array([zstart])  
    # pdb.set_trace()
    edges_src = mlab.pipeline.scalar_scatter(
        xv.ravel(), yv.ravel(), zv.ravel(), edge_c.ravel())

    edges_src.mlab_source.dataset.point_data.add_array(
        np.abs(np.array(edge_r).ravel()))
    edges_src.mlab_source.dataset.point_data.get_array(1).name = 'radius'

    cell_array.set_cells(edgecount, lines)
    
    edges_src.mlab_source.dataset.lines = cell_array

    edges_src.mlab_source.update()

    if tube_radius is not None:
        #edges_src.mlab_source.dataset.point_data.vectors = np.abs(
        #        np.array([edge_r,edge_r,edge_r]).T)
        edges_src.mlab_source.dataset.point_data.add_array(
                np.abs(np.array(edge_r).ravel()))
        edges_src.mlab_source.dataset.point_data.get_array(1).name = 'radius'
        tubes = mlab.pipeline.tube(mlab. pipeline.set_active_attribute(edges_src,
                                    point_scalars='radius'),
                                     tube_radius=tube_radius,)
        
        #tubes.filter.radius=edge_r
        tubes.filter.vary_radius = 'vary_radius_by_scalar'
    else:
        tubes = edges_src
    tube_surf = mlab.pipeline.surface(mlab.pipeline. set_active_attribute(
        tubes, point_scalars='scalars'), colormap=edge_colormap, color=edge_color, **kwargs)
    tubes.children[0].point_scalars_name='scalars'
    if edge_color_keyword is None:
        tube_surf.actor.mapper.scalar_visibility = False
    else:
        tube_surf.actor.mapper.scalar_visibility = True
    #pdb.set_trace()    
    if not node_size:
        return tube_surf, None

    # ------------------------------------------------------------------------
    # Plot the nodes
    if node_size is not None:
        nodes = G.nodes()
        L = len(G.nodes())
        
        if node_color_keyword is None and scale_node_keyword is None:
            
            nodes = G.nodes()
            L = len(G.nodes())

           
            s = np.ones(len(xn))
            # pdb.set_trace()
            pts = mlab.points3d(xn, yn, zn, s,
                                scale_factor=node_size,
                                # colormap=node_colormap,
                                resolution=16,
                                color=node_color)

        else:
         
            if  scale_node_keyword is not None:
               if scale_node_keyword == 'degree':
                    dic = G.degree(nodes)
                    s = []
                    for n in nodes:
                            s.append(dic[n])
               else:
                    s= nx.get_node_attributes(G,scale_node_keyword).values() 
            else:
                s = np.ones(len(xn)) * node_size   
                      
            if node_color_keyword is not None:
                node_color_scalar = nx.get_node_attributes(
                    G, node_color_keyword ).values()
            else:
                node_color_scalar = len(xn) * [1]

            pts = mlab.quiver3d(np.array(xn).ravel(), np.array(
                yn).ravel(), np.array(zn).ravel(),s,s,s, scalars=np.array(node_color_scalar).ravel(), mode='sphere',scale_factor=node_size,resolution=16)
            
            pts.glyph.color_mode = 'color_by_scalar'
            #pts.glyph.scale_mode = 'scale_by_vector'
    print 'done'
    fig.scene.disable_render = disable_render
    return tube_surf, pts


