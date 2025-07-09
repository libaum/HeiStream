#ifndef MAP_STORAGE_H
#define MAP_STORAGE_H

#include "pq_storage_interface.h"
#include "pq_item.h"
#include <unordered_map>

// Map-basierte Implementierung
template <typename Key = LongNodeID, typename Value = PQItem>
class MapStorage : public PQStorageBase<Key, Value> {
    std::unordered_map<Key, Value> data;

public:
    void insert(const Key& key, const Value& value) override {
        data[key] = value;
    }

    void emplace(const Key& key, Value&& value) override {
        data.emplace(key, std::move(value));
    }

    std::optional<std::reference_wrapper<Value>> get(const Key& key) override {
        auto it = data.find(key);
        if (it != data.end()) return std::ref(it->second);
        return std::nullopt;
    }

    std::optional<std::reference_wrapper<const Value>> get(const Key& key) const override {
        auto it = data.find(key);
        if (it != data.end()) return std::cref(it->second);
        return std::nullopt;
    }

    void erase(const Key& key) override {
        data.erase(key);
    }

    bool contains(const Key& key) const override {
        return data.find(key) != data.end();
    }

    Value& operator[](const Key& key) override {
        return data[key];
    }

    void reserve(size_t size) override {
        data.reserve(size);
    }

    void set_max_load_factor(float factor) override {
        data.max_load_factor(factor);
    }
};

#endif
