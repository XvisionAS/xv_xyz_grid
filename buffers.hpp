#ifndef __BUFFERS_H__
#define __BUFFERS_H__

#include <map>
#include <vector>

template <typename T>
struct tBuffer {
	std::vector<T>	mData;
};


template <typename T>
struct tBuffers {
	static tBuffers& singleton() {
		static tBuffers singleton;
		return singleton;
	}

	void setBufferLimit(const size_t memoryLimit) {
		mMemoryLimit = memoryLimit;
		mBuffers.clear();
	}

	tBuffer<T>& get(const std::string& fileName) {
		auto it = mBuffers.find(fileName);	
		if (it != mBuffers.end()) {
			return it->second;
		}
		else {
			discardIfNeeded();
			auto& ret = mBuffers[fileName];
			fromFile(fileName, ret.mData);
			return ret;
		}		
	}

	void discardIfNeeded() {
		size_t totalMemory = 0;
		for (auto m = mBuffers.begin(); m != mBuffers.end(); ) {
			totalMemory += m->second.mData.size() * sizeof(T);
			if (totalMemory > mMemoryLimit) {
				toFile(m->first, m->second.mData);
				mBuffers.erase(m++);
			} else {
				++m;
			}
		}
	}

	void discardAll() {
		for (auto m = mBuffers.begin(); m != mBuffers.end(); ++m) {
			toFile(m->first, m->second.mData);
		}
		mBuffers.clear();
	}

private:
	tBuffers() : mMemoryLimit( 100 * 1000 * 1000 ) {}
	tBuffers(const tBuffers&) {}

	size_t										mMemoryLimit;
	std::map< std::string, tBuffer<T> >	mBuffers;

	void toFile(const std::string& fileName, const std::vector<T>& data) {
		std::ofstream out;
		out.open(fileName.c_str(), std::ios::out | std::ios::binary);
		int size = data.size();
		out.write((char*)&size, sizeof(size));
		out.write((char*)&data[0], size * sizeof(T));
	}

	void fromFile(const std::string& fileName, std::vector<T>& data) {
		std::ifstream out;
		out.open(fileName.c_str(), std::ios::in | std::ios::binary);
		if (out.is_open()) {
			int size;
			out.read((char*)&size, sizeof(size));
			data.resize(size);
			out.read((char*)&data[0], size * sizeof(T));
		}
	}
};


#endif
