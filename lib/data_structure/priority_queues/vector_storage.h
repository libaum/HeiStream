#ifndef VECTOR_STORAGE_H
#define VECTOR_STORAGE_H

#include "pq_storage_interface.h"
#include "pq_item.h"
#include <vector>

// Vector-basierte Implementierung
template <typename Key = LongNodeID, typename Value = PQItem>
class VectorStorage : public PQStorageBase<Key, Value> {
    std::vector<Value*> data;

public:
    VectorStorage() = default;

    VectorStorage(size_t initial_size) : data(initial_size, nullptr) {}

    ~VectorStorage() {
        for (auto* ptr : data) {
            delete ptr;
        }
    }

    void insert(const Key& key, const Value& value) override {
        if (data[key - 1] != nullptr) {
            delete data[key - 1];
        }
        data[key - 1] = new Value(value);
    }

    void emplace(const Key& key, Value&& value) override {
        if (data[key - 1] != nullptr) {
            delete data[key - 1];
        }
        data[key - 1] = new Value(std::move(value));
    }

    std::optional<std::reference_wrapper<Value>> get(const Key& key) override {
        if (data[key - 1] != nullptr) {
            return std::ref(*data[key - 1]);
        }
        return std::nullopt;
    }

    std::optional<std::reference_wrapper<const Value>> get(const Key& key) const override {
        if (data[key - 1] != nullptr) {
            return std::cref(*data[key - 1]);
        }
        return std::nullopt;
    }

    void erase(const Key& key) override {
        if (data[key - 1] != nullptr) {
            delete data[key - 1];
            data[key - 1] = nullptr;
        }
    }

    bool contains(const Key& key) const override {
        return data[key - 1] != nullptr;
    }

    Value& operator[](const Key& key) override {
        if (data[key - 1] == nullptr) {
            data[key - 1] = new Value{};
        }
        return *data[key - 1];
    }

    void reserve(size_t size) override {
        // data.reserve(size);
    }

    void set_max_load_factor(float factor) override {
        // Vector hat keinen load factor, ignoriere
    }
};

#endif
