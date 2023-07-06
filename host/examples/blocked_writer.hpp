#ifndef __BLOCKED_WRITER
#define __BLOCKED_WRITER

#include <iostream>
#include <memory>

static constexpr int num_writes = 1000; 

class blocked_writer {
private:
    std::unique_ptr<std::ofstream> of;
    std::string file_name;
    int num = 0, index = 0;
public:
    blocked_writer(const std::string& file_name): file_name(file_name) {
    }

    void write(const char* buff, size_t n) {
        if (num == 0) {
            of = std::make_unique<std::ofstream>(file_name + "_" + std::to_string(index++), std::ofstream::binary);
        }
        of->write(buff, n);
        if (++num == num_writes) {
            num = 0;
            of->close();
        }
    }
};

#endif