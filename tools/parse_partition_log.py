#!/usr/bin/env python3
"""
Utility to inspect BuffCut FlatBuffer log files produced with --write_log.

The script only depends on the Python flatbuffers runtime:
    pip install flatbuffers

Usage:
    python tools/parse_partition_log.py path/to/log.bin
"""

import argparse
from pathlib import Path
from typing import Optional

try:
    import flatbuffers
    from flatbuffers.table import Table
except ImportError as exc:  # pragma: no cover - guidance for users
    raise SystemExit(
        "flatbuffers package missing. Install it via 'pip install flatbuffers'."
    ) from exc


class GraphMetadata:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def Filename(self) -> str:
        o = self._tab.Offset(4)
        if o:
            return self._tab.String(o + self._tab.Pos).decode("utf-8")
        return ""

    def NumNodes(self) -> int:
        o = self._tab.Offset(6)
        return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos) if o else 0

    def NumEdges(self) -> int:
        o = self._tab.Offset(8)
        return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos) if o else 0


class PartitionConfiguration:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def K(self) -> int:
        o = self._tab.Offset(4)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def Seed(self) -> int:
        o = self._tab.Offset(6)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def BatchSize(self) -> int:
        o = self._tab.Offset(8)
        return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos) if o else 0

    def ModelMode(self) -> int:
        o = self._tab.Offset(10)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def Alpha(self) -> int:
        o = self._tab.Offset(12)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def MaxBufferSize(self) -> int:
        o = self._tab.Offset(14)
        return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos) if o else 0

    def BbRatio(self) -> int:
        o = self._tab.Offset(16)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def BufferScoreType(self) -> int:
        o = self._tab.Offset(18)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def DMax(self) -> int:
        o = self._tab.Offset(20)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def HaaBeta(self) -> float:
        o = self._tab.Offset(22)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def HaaTheta(self) -> float:
        o = self._tab.Offset(24)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def CbsTheta(self) -> float:
        o = self._tab.Offset(26)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def BufferNeighborWeight(self) -> float:
        o = self._tab.Offset(28)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def BatchFrontierWeight(self) -> float:
        o = self._tab.Offset(30)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def NumStreamsPasses(self) -> int:
        o = self._tab.Offset(32)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def RestreamNumber(self) -> int:
        o = self._tab.Offset(34)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def RestreamIncludeHighDegree(self) -> bool:
        o = self._tab.Offset(36)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False

    def RestreamVcycle(self) -> bool:
        o = self._tab.Offset(38)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False

    def GhostNeighborsEnabled(self) -> bool:
        o = self._tab.Offset(40)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False

    def GhostWeight(self) -> float:
        o = self._tab.Offset(42)
        return self._tab.Get(flatbuffers.number_types.Float32Flags, o + self._tab.Pos) if o else 0.0

    def SepBatchMarker(self) -> bool:
        o = self._tab.Offset(44)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False

    def BatchExtractionStrategy(self) -> int:
        o = self._tab.Offset(46)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def MaxActiveBatches(self) -> int:
        o = self._tab.Offset(48)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def MaxInputQueueSize(self) -> int:
        o = self._tab.Offset(50)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def AltThreadQueue(self) -> bool:
        o = self._tab.Offset(52)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False

    def CollectLocalityMetrics(self) -> bool:
        o = self._tab.Offset(54)
        return bool(self._tab.Get(flatbuffers.number_types.Uint8Flags, o + self._tab.Pos)) if o else False


class RunTime:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def IoTime(self) -> float:
        o = self._tab.Offset(4)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def PartitionTime(self) -> float:
        o = self._tab.Offset(6)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def ModelConstructionTime(self) -> float:
        o = self._tab.Offset(8)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def MappingTime(self) -> float:
        o = self._tab.Offset(10)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def TotalTime(self) -> float:
        o = self._tab.Offset(12)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0


class MemoryConsumption:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def MaxRSS(self) -> int:
        o = self._tab.Offset(4)
        return self._tab.Get(flatbuffers.number_types.Int64Flags, o + self._tab.Pos) if o else 0


class PartitionMetrics:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def EdgeCut(self) -> int:
        o = self._tab.Offset(4)
        return self._tab.Get(flatbuffers.number_types.Int32Flags, o + self._tab.Pos) if o else 0

    def VertexCut(self) -> int:
        o = self._tab.Offset(6)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def Replicas(self) -> int:
        o = self._tab.Offset(8)
        return self._tab.Get(flatbuffers.number_types.Uint32Flags, o + self._tab.Pos) if o else 0

    def ReplicationFactor(self) -> float:
        o = self._tab.Offset(10)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def Balance(self) -> float:
        o = self._tab.Offset(12)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def InternalEdgeRatio(self) -> float:
        o = self._tab.Offset(14)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def NeighborCoverage(self) -> float:
        o = self._tab.Offset(16)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0

    def Conductance(self) -> float:
        o = self._tab.Offset(18)
        return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos) if o else 0.0


class PartitionLog:
    __slots__ = ("_tab",)

    def Init(self, buf: bytearray, pos: int) -> None:
        self._tab = Table(buf, pos)

    def GraphMetadata(self) -> Optional[GraphMetadata]:
        o = self._tab.Offset(4)
        if o:
            obj = GraphMetadata()
            obj.Init(self._tab.Bytes, self._tab.Indirect(o + self._tab.Pos))
            return obj
        return None

    def PartitionConfiguration(self) -> Optional[PartitionConfiguration]:
        o = self._tab.Offset(6)
        if o:
            obj = PartitionConfiguration()
            obj.Init(self._tab.Bytes, self._tab.Indirect(o + self._tab.Pos))
            return obj
        return None

    def RunTime(self) -> Optional[RunTime]:
        o = self._tab.Offset(8)
        if o:
            obj = RunTime()
            obj.Init(self._tab.Bytes, self._tab.Indirect(o + self._tab.Pos))
            return obj
        return None

    def MemoryConsumption(self) -> Optional[MemoryConsumption]:
        o = self._tab.Offset(10)
        if o:
            obj = MemoryConsumption()
            obj.Init(self._tab.Bytes, self._tab.Indirect(o + self._tab.Pos))
            return obj
        return None

    def Metrics(self) -> Optional[PartitionMetrics]:
        o = self._tab.Offset(12)
        if o:
            obj = PartitionMetrics()
            obj.Init(self._tab.Bytes, self._tab.Indirect(o + self._tab.Pos))
            return obj
        return None


def get_root(buf: bytearray) -> PartitionLog:
    obj = PartitionLog()
    obj.Init(buf, flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, 0))
    return obj


def print_partition_log(log: PartitionLog) -> None:
    metadata = log.GraphMetadata()
    config = log.PartitionConfiguration()
    runtime = log.RunTime()
    memory = log.MemoryConsumption()
    metrics = log.Metrics()

    print("== Graph Metadata ==")
    if metadata:
        print(f"Graph:      {metadata.Filename()}")
        print(f"Nodes (n):  {metadata.NumNodes()}")
        print(f"Edges (m):  {metadata.NumEdges()}")

    print("\n== Partition Configuration ==")
    if config:
        print(f"k blocks:   {config.K()}")
        print(f"Seed:       {config.Seed()}")
        print(f"Batch size: {config.BatchSize()}")
        model_mode = config.ModelMode()
        alpha = config.Alpha()
        print(f"Model mode: {model_mode}")
        print(f"Alpha:      {alpha}")
        print(f"Max buffer: {config.MaxBufferSize()}")
        print(f"BB ratio:   {config.BbRatio()}")
        print(f"Buffer score type: {config.BufferScoreType()}")
        print(f"d_max:      {config.DMax()}")
        print(f"HAA beta:   {config.HaaBeta():.3f}")
        print(f"HAA theta:  {config.HaaTheta():.3f}")
        print(f"CBS theta:  {config.CbsTheta():.3f}")
        print(f"buffer_neighbor_w: {config.BufferNeighborWeight():.3f}")
        print(f"batch_frontier_w:  {config.BatchFrontierWeight():.3f}")
        print(f"Streams passes:    {config.NumStreamsPasses()}")
        print(f"Restream number:   {config.RestreamNumber()}")
        print(f"Restream include high degree: {config.RestreamIncludeHighDegree()}")
        print(f"Restream vcycle:   {config.RestreamVcycle()}")
        print(f"Ghost neighbors:   {config.GhostNeighborsEnabled()} (weight={config.GhostWeight():.3f})")
        print(f"Separate batch markers: {config.SepBatchMarker()}")
        print(f"Batch extraction strategy: {config.BatchExtractionStrategy()}")
        print(f"Max active batches: {config.MaxActiveBatches()}")
        print(f"Max input queue:    {config.MaxInputQueueSize()}")
        print(f"Alt thread queue:   {config.AltThreadQueue()}")
        print(f"Collect locality metrics: {config.CollectLocalityMetrics()}")

    print("\n== Runtime ==")
    if runtime:
        print(f"IO time:                {runtime.IoTime():.6f} s")
        print(f"Model construction:     {runtime.ModelConstructionTime():.6f} s")
        print(f"Mapping time:           {runtime.MappingTime():.6f} s")
        print(f"Partition time:         {runtime.PartitionTime():.6f} s")
        print(f"Total time:             {runtime.TotalTime():.6f} s")

    print("\n== Memory ==")
    if memory:
        print(f"Peak RSS: {memory.MaxRSS()} KB")

    print("\n== Partition Metrics ==")
    if metrics:
        print(f"Edge cut:           {metrics.EdgeCut()}")
        print(f"Vertex cut:         {metrics.VertexCut()}")
        print(f"Replicas:           {metrics.Replicas()}")
        print(f"Replication factor: {metrics.ReplicationFactor():.6f}")
        print(f"Balance:            {metrics.Balance():.6f}")
        print(f"Avg internal edges: {metrics.InternalEdgeRatio():.6f}")
        print(f"Avg neighbor cover: {metrics.NeighborCoverage():.6f}")
        print(f"Avg conductance:    {metrics.Conductance():.6f}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect a BuffCut FlatBuffer log.")
    parser.add_argument("bin_path", type=Path, help="Path to the <graph>_k_batch_buffer.bin file")
    args = parser.parse_args()

    data = bytearray(args.bin_path.read_bytes())
    log = get_root(data)
    print_partition_log(log)


if __name__ == "__main__":
    main()
