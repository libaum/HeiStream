
# HeiStream Enhanced Version

## Overview



This repository contains an enhanced version of the HeiStream algorithm, a buffered streaming algorithm designed to solve the graph partitioning problem. Our enhancements aim to improve the robustness and performance of HeiStream, particularly by addressing its dependency on node ordering during streaming.
For more general information about HeiStream have a look at the original repository.

### Key Enhancements

1.  **Global Approach**: Implements a delayed batch insertion queue for nodes with low locality, reducing the impact of suboptimal node ordering (in *graph_io_stream::loadLinesFromStreamToBinary*).

2.  **Local Approach**: Introduces a priority queue for initial partitioning to enhance the partition order (in *init_fennel::fennel*).



These modifications aim to improve partitioning quality and efficiency, especially under random or suboptimal node orderings.


## Running the Enhanced Version of HeiStream


To partition a graph in METIS format using the basic model of HeiStream, run


```shell
./heistream <graph  filename> --k=<number  of  blocks> --stream_buffer=<nodes  per  bufer>
```
For our global approach, we added the optional parameters (both default to 0):
- `--max_delayed_nodes=<maximum  number  of  delayed  nodes>`

- `--threshold_delay=<delay  threshold>`


To activate the local approach use the parameter:
- `--local_pq`


### Examples

 To execute the only the global approach, this would be a valid command:

```shell
./heistream  graphs/in-2004.graph  --k=2  --max_delayed_nodes=131072  --threshold_delay=0.1
```

  To execute only the local approach, this would be a possible command:

```shell
./heistream  graphs/in-2004.graph  --k=2 --local_pq"
```

