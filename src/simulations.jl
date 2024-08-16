# Simulation of branching process 
# This is only used to make the tree in Figure 1 and none of the scientific conclusions rely on it. 


"""
    Cell{T}

A mutable struct representing a cell in a tree.

## Fields

- `label::T`: Data associated with the cell.
- `left::Union{Cell{T}, Nothing}`: Left child cell.
- `right::Union{Cell{T}, Nothing}`: Right child cell.

"""
mutable struct Cell{T}
    label::T
    left::Union{Cell{T}, Nothing}
    right::Union{Cell{T}, Nothing}
end

"""
    create_cell(label, left=nothing, right=nothing)

Create a new cell with the given label and optional left and right children.

## Arguments

- `label`: Data to be stored in the cell.
- `left=nothing`: Left child cell (default is `nothing`).
- `right=nothing`: Right child cell (default is `nothing`).

## Returns

- `Cell`: A new `Cell` object with the specified label and children.

"""
function create_cell(label, left=nothing, right=nothing)
    return Cell(label, left, right)
end


function makelineage(generator,init,n,params,labels)
    cells = [generator(init, params)]
    for _ in 2:n
        push!(cells,generator(cells[end],params))
    end
    
    df = DataFrame(hcat([c for c in cells]...)',labels)
    df[:,:n] = cumsum(1:length(df[:,1]))
    df
end





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# POPULATION GROWTH SIMULATIONS
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
"""
    grow_tree!(cell, terminate, params, generator)

Recursively grow a tree from a single root cell.

## Arguments

- `cell`: The root cell from which to grow the tree.
- `terminate`: A function to determine when to terminate growth.
- `params`: Parameters used in the growth process.
- `generator`: A function to generate new cell labels.

"""
function grow_tree!(cell, terminate, params, generator)
    L = generator(cell.label, params)
    R = generator(cell.label, params)
    if terminate(cell) == false
        cell.left = create_cell(L)
        cell.right = create_cell(R)
        grow_tree!(cell.left, terminate, params, generator)
        grow_tree!(cell.right, terminate, params, generator)
    end
end


"""
    grow_forest!(nodes, terminate, params, dist)

Grow a forest of trees from a collection of root nodes.

## Arguments

- `nodes`: Collection of root nodes for trees.
- `terminate`: A function to determine when to terminate growth.
- `params`: Parameters used in the growth process.
- `dist`: Distribution used for generating new cell labels.

"""
function grow_forest!(nodes, terminate, params, dist)
    for node in nodes
        grow_tree!(node, terminate, params, dist)
    end
end

"""
    find_max(node::Cell)

Recursively find the maximum value in a tree of cells.

## Arguments

- `node::Cell`: The root node of the tree.

## Returns

- `T`: The maximum value in the tree.

"""
function find_max(node::Cell)
    if node.left != nothing
        left_max = find_max(node.left)
        right_max = find_max(node.right)
        current_max = max(node.label, max(left_max, right_max))
        return current_max
    else
        return node.label
    end
end

"""
    get_leaf_nodes(node)

Get a list of leaf nodes in a tree.

## Arguments

- `node`: The root node of the tree.

## Returns

- `Array`: An array of leaf nodes

"""
function get_leaf_nodes(node)
    if (node.left != nothing) && (node.right != nothing)
        L = get_leaf_nodes(node.left)
        R = get_leaf_nodes(node.right)
        return vcat(L, R)
    else
        return [node]
    end
end

"""
    get_lineage(node)

Get the lineage of a node in the tree.

## Arguments

- `node`: The node for which to retrieve the lineage.

## Returns

- `Array`: An array of node labels representing the lineage.

"""
function get_lineage(node)
    lineage = [node.label]
    while (node.left != nothing) && (node.right != nothing)
        r = rand()
        if r < 0.5
            node = node.left
        else
            node = node.right
        end
        push!(lineage, node.label)
    end
    return vcat(lineage...)
end

"""
    get_node(root, path)

Get a node from the tree given a path.

## Arguments

- `root`: The root node of the tree.
- `path`: An array of 0s and 1s representing the path to the desired node.

## Returns

- `T`: The label of the node at the specified path.

"""
function get_node(root, path)
    node = root
    for p in path[0:-1]
        if p == 0
            node = node.left
        elseif p == 1
            node = node.right
        end
    end
    return node.label
end

