/******************************************************************************
 * flat_buffer_writer.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adilchhabra7@gmail.com>
 *****************************************************************************/

 #ifndef KAHIP_FLATBUFFERWRITER_H
 #define KAHIP_FLATBUFFERWRITER_H

 #include <fstream>
 #include <iostream>
 #include <vector>
 #include "PartitionInfo_generated.h"
 #include "partition/partition_config.h"
 #include <iomanip>

 const bool SUPPRESS_OUTPUT = true;

 class FlatBufferWriter {
 private:
     double io_time_;
     double model_construction_time_;
     double mapping_time_;
     double partition_time_;
     double total_time_;
     long maxRSS_;
     EdgeWeight edge_cut_;
     NodeID vertex_cut_;
     NodeID replicas_;
     double replication_factor_ ;
     double balance_;

 public:
     FlatBufferWriter()
             : io_time_(0.0), model_construction_time_(0.0), mapping_time_(0.0),
               partition_time_(0.0), total_time_(0.0), maxRSS_(0), edge_cut_(-1),
               vertex_cut_(-1), replicas_(-1), replication_factor_(0.0), balance_ (0.0){}

     void updateResourceConsumption(double &io_time, double &model_construction_time,
                                     double &mapping_time, double &partition_time, double &total_time, long &maxRSS) {
         io_time_ = io_time;
         model_construction_time_ = model_construction_time;
         mapping_time_ = mapping_time;
         partition_time_ = partition_time;
         total_time_ = total_time;
         maxRSS_ = maxRSS;
     }

     void updateVertexPartitionResults(EdgeWeight & edge_cut, double balance) {
         edge_cut_ = edge_cut;
         balance_ = balance;
     }

     void updateEdgePartitionResults(NodeID & vertex_cut, NodeID & replicas, double replication_factor, double balance) {
         vertex_cut_ = vertex_cut;
         replicas_ = replicas;
         replication_factor_ = replication_factor;
         balance_ = balance;
     }

     // Function to extract the base filename without path and extension
     static std::string extractBaseFilename(const std::string& fullPath) {
         size_t lastSlash = fullPath.find_last_of('/');
         size_t lastDot = fullPath.find_last_of('.');

         if (lastSlash != std::string::npos) {
             // Found a slash, extract the substring after the last slash
             return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
         } else {
             // No slash found, just extract the substring before the last dot
             return fullPath.substr(0, lastDot);
         }
     }

     void write(const std::string& graph_filename, const PartitionConfig &partition_config ) const {
         flatbuffers::FlatBufferBuilder builder(1024);
         std::stringstream filename;
         std::string baseFilename = extractBaseFilename(graph_filename);
         auto filenameOffset = builder.CreateString(baseFilename);

         auto metadata = PartitionInfo::CreateGraphMetadata(
            builder, filenameOffset, partition_config.total_nodes, partition_config.total_edges);

        auto configdata = PartitionInfo::CreatePartitionConfiguration(
            builder, partition_config.k, partition_config.seed, partition_config.stream_buffer_len,
            partition_config.minimal_mode ? -1 : 0,
            partition_config.batch_alpha ? 0 : 0);

         if(partition_config.edge_partition) {
             if (partition_config.minimal_mode) std::cout << "Mode: Minimal" << std::endl;
             if (partition_config.batch_alpha) std::cout << "Alpha: Batch Alpha" << std::endl;
         }

         auto runtimedata = PartitionInfo::CreateRunTime(
            builder, io_time_, model_construction_time_, mapping_time_,
            partition_time_, total_time_);

        EdgeWeight edge_cut_val = partition_config.edge_partition ? 0 : edge_cut_;
        NodeID vertex_cut_val = partition_config.edge_partition ? vertex_cut_ : 0;
        NodeID replicas_val = partition_config.edge_partition ? replicas_ : 0;
        double replication_factor_val = partition_config.edge_partition ? replication_factor_ : 0;

        auto partition_metrics = PartitionInfo::CreatePartitionMetrics(
            builder, edge_cut_val, vertex_cut_val, replicas_val,
            replication_factor_val, balance_);

         // Create MemoryConsumption
         auto memory_consumption = PartitionInfo::CreateMemoryConsumption(builder, maxRSS_);

         // Create PartitionLog
         auto partition_log = PartitionInfo::CreatePartitionLog(
            builder, metadata, configdata, runtimedata, memory_consumption, partition_metrics);

        builder.Finish(partition_log);

         //Step 4: Write to File
         const uint8_t* bufferPointer = builder.GetBufferPointer();
         int bufferSize = builder.GetSize();

         std::string outputFileNameStream;
         outputFileNameStream = baseFilename + "_" + std::to_string(partition_config.k) + "_" + std::to_string(partition_config.stream_buffer_len) + "_" + std::to_string(partition_config.max_pq_size) + ".bin";
         const char* outputFileName = outputFileNameStream.c_str();
         if(partition_config.write_log) {
             FILE *file = fopen(outputFileName, "wb");
             fwrite(bufferPointer, 1, bufferSize, file);
             fclose(file);
         }

         if (!SUPPRESS_OUTPUT){
             std::cout << "Blocks = k: " << partition_config.k << std::endl;
             std::cout << "Seed: " << partition_config.seed << std::endl;
             std::cout << "Stream Buffer: " << partition_config.stream_buffer_len << std::endl;
             std::cout << "Graph: " << baseFilename << std::endl;
             std::cout << "Nodes (n): " << partition_config.total_nodes << std::endl;
             std::cout << "Edges (m): " << partition_config.total_edges << std::endl;
             std::cout << "IO time: " << io_time_ << std::endl;
             std::cout << "Model Construction time: " << model_construction_time_ << std::endl;
             std::cout << "Mapping time: " << mapping_time_ << std::endl;
             std::cout << "Partition time: " << partition_time_ << std::endl;
             std::cout << "Total time: " << total_time_ << std::endl;
             if(!partition_config.edge_partition) {
                 std::cout << "Edge Cut: " << edge_cut_ << std::endl;
             }
             std::cout << "Balance: " << balance_ << std::endl;
             if (maxRSS_ != -1) {
                 std::cout << "Maximum Resident Set Size (KB): " << maxRSS_ << std::endl;
             }
         }
     }
 };


 #endif //KAHIP_FLATBUFFERWRITER_H
