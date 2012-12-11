//
// Kmernator/src/LRUCache.h
//

#ifndef LRU_CACHE_H
#define LRU_CACHE_H

#include <mct/hash-map.hpp>

template<typename KEY, typename VALUE>
class LRUCache {
public:
	typedef KEY KeyType;
	typedef VALUE ValueType;
	typedef typename mct::linked_hash_map<KEY, VALUE> Impl;
	typedef typename Impl::iterator ImplIterator;

	LRUCache(const unsigned long _maxSize) : maxSize(_maxSize){

	}
	bool fetch(const KeyType &key, ValueType &val) {
		ImplIterator entry = cache.find(key);
		if (entry != cache.end()) {
			val = entry->second;
			cache.relink(cache.end(), entry);
			return true;
		} else {
  			return false;
		}
	}
	void insert(const KeyType &key, const ValueType &val) {
		if (cache.size() >= maxSize)
			cache.pop_front();
		cache.insert( std::make_pair(key, val) );
	}
private:
	Impl cache;
	unsigned long maxSize;
};

#endif // LRRU_CACHE_H
