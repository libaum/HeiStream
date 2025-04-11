// Im Header-Bereich oder in einer neuen Datei mlp_thread_manager.h
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>
#include <memory>
#include "partition/partition_config.h"
#include "graph_io_stream.h"
#include "partition/graph_partitioner.h"

void config_multibfs_initial_partitioning(PartitionConfig &partition_config);

void config_multibfs_initial_partitioning(PartitionConfig &partition_config) {
    if (partition_config.initial_part_multi_bfs && partition_config.curr_batch >= 2) {
        partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
    }
}

// A function to do multi-level partitioning the nodes in the batch (input)
void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> *&batch_nodes) {
// void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<LongNodeID> *&input_idxs, Buffer &buffer) {
    // Initialize the partition configuration
    graph_access G = graph_access();
    quality_metrics qm;
    balance_configuration bc;

    // ***************************** build model ***************************************
    G.set_partition_count(partition_config.k);
    partition_config.local_to_global_map = new std::vector<NodeID>(partition_config.nmbNodes);
    graph_io_stream::createModel(partition_config, G, batch_nodes);
    graph_io_stream::countAssignedNodes(partition_config);
    graph_io_stream::prescribeBufferInbalance(partition_config);
    bool already_fully_partitioned = (partition_config.restream_vcycle && partition_config.restream_number);
    bc.configurate_balance(partition_config, G, already_fully_partitioned || !partition_config.stream_initial_bisections);
    config_multibfs_initial_partitioning(partition_config);

    // ***************************** perform partitioning ***************************************
    graph_partitioner partitioner;
    partitioner.perform_partitioning(partition_config, G);

    // ***************************** permanent assignment ***************************************
    graph_io_stream::generalizeStreamPartition(partition_config, G);

    delete partition_config.local_to_global_map;
}

// Thread-Management für MLP
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

class MLPThreadManager {
private:
    std::thread worker;
    std::mutex mtx;
    std::condition_variable task_cv;
    std::condition_variable completion_cv;

    std::atomic<bool> is_running{false};
    bool has_task = false;
    bool shutdown = false;

    PartitionConfig* config_ptr = nullptr;
    std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>** batch_ptr = nullptr;

    void worker_loop() {
        while (!shutdown) {
            // Auf Aufgabe warten
            {
                std::unique_lock<std::mutex> lock(mtx);
                task_cv.wait(lock, [this]{ return has_task || shutdown; });

                if (shutdown) break;

                is_running = true;
            }

            // Ausführen der MLP-Funktion
            if (config_ptr && batch_ptr) {
                ::perform_mlp_on_batch(*config_ptr, *batch_ptr);
            }

            // Signalisiere Fertigstellung
            {
                std::lock_guard<std::mutex> lock(mtx);
                has_task = false;
                is_running = false;
            }
            completion_cv.notify_all();
        }
    }

public:
    MLPThreadManager() {
        worker = std::thread(&MLPThreadManager::worker_loop, this);
    }

    ~MLPThreadManager() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            shutdown = true;
        }
        task_cv.notify_one();
        if (worker.joinable()) {
            worker.join();
        }
    }

    // Führt perform_mlp_on_batch asynchron aus (wartet wenn nötig)
    void execute(PartitionConfig& config,
                std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>*& batch) {
        std::unique_lock<std::mutex> lock(mtx);
        // Warten wenn bereits eine Aufgabe läuft
        completion_cv.wait(lock, [this]{ return !has_task; });

        // Parameter setzen und Thread starten
        config_ptr = &config;
        batch_ptr = &batch;
        has_task = true;

        lock.unlock();
        task_cv.notify_one();
    }

    // Warten bis die aktuelle MLP-Berechnung abgeschlossen ist
    void wait_completion() {
        std::unique_lock<std::mutex> lock(mtx);
        completion_cv.wait(lock, [this]{ return !is_running && !has_task; });
    }

    bool is_busy() const {
        return is_running.load();
    }
};

