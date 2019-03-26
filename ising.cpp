// g++ --std=c++14 -I. -o ising -O3 ising.cpp # or clang++
#include <matrix.h>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <cwchar>
#include <functional>
using namespace std;

class World: public Matrix
{
    double beta_;
    mt19937 generator_;
    mt19937 generator2_;
    uniform_int_distribution<int> row_picker_;
    uniform_int_distribution<int> col_picker_;
    uniform_real_distribution<double> dist_;
    function<double()> rnd;
    function<int()> rnd_row;
    function<int()> rnd_col;
    
public:
    using Matrix::Matrix; // c++11: matrix constructors
    World(uint32_t rows, uint32_t cols, double temp=1.0, int seed=0)
            : Matrix(rows, cols), beta_(1./temp),
              row_picker_(0, rows - 1), col_picker_(0, cols - 1)
    {
        generator_.seed(seed);
        generator2_.seed(seed + 1);
        rnd = bind(dist_, generator_);
        rnd_row = bind(row_picker_, generator_);
        rnd_col = bind(col_picker_, generator2_);
    }

    void init(double fraction, int seed=0)
    {
        bernoulli_distribution dist(fraction);
        generator_.seed(seed);
        auto rnd_init = bind(dist, generator_);
        generate(begin(data()), end(data()), [&rnd_init]() { return rnd_init() ? 1 : -1; });
    }

    void set_temp(double temp)
    {
        beta_ = 1./temp;
    }

    int neighbour_sum(int row, int col) const
    {
        auto R = getRows();
        auto C = getCols();
        return get((row + 1) % R, col) +
            get((row - 1 + R) % R, col) +
            get(row, (col + 1) % C) +
            get(row, (col - 1 + C) % C);
    }
    
    bool update(int n=1) 
    {
        bool changed = false;
        
        for (int i = 0; i < n; i++)
        {
            int row = rnd_row();
            int col = rnd_col();
            auto val = get(row, col);
            double delta_E = 2 * val * neighbour_sum(row, col);
            // printf("(%d,%d), ΔE = %.1f: ", row, col, delta_E); ///
            // accept or reject?
            if (rnd() < exp(-beta_ * delta_E))
            {
                set(row, col, -val);
                changed = true;
                // printf("accepted\n"); ///
            }
            else
            {
                // printf("rejected\n"); ///
            }
        }
        return changed;
    }
    
    void print() const
    {
        uint32_t R = getRows();
        uint32_t C = getCols();
        for (uint32_t r = 0; r < R; r++)
        {
            for (uint32_t c = 0; c < C; c++)
            {
                putchar(get(r,c) == 1 ? 'O' : ' ');
                //putwchar(get(r,c) == 1 ? u'↑' : u'↓');
            }
            if (r != R - 1)
            {
                putchar('\n');
            }
        }
        cout << flush;
    }
};
    

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <sys/ioctl.h>
#include <unistd.h>
#include <fcntl.h>
#include <linux/fb.h>
#include <sys/mman.h>

int main_txt(int generations, int steps_per_generation,
             int delay, double fraction, int seed, double temp)
{
    struct winsize size;
    ioctl(STDOUT_FILENO,TIOCGWINSZ,&size);

    World m(size.ws_row, size.ws_col, temp, seed);
    // World m(5,13, temp);
    // printf("Size: %d x %d\n", size.ws_row, size.ws_col);
    // exit(0);
    m.init(fraction);
    m.print();

    for (int i = 0; i < generations; i++)
    {
        this_thread::sleep_for(std::chrono::milliseconds(delay));
        m.update(steps_per_generation);
        m.print();
        // printf("\n");printf("-------------\n"); ///
    }
    return 0;
}

// fbfd is the open (R/W) file descriptor of the framebuffer
int main_fb(int fbfd, int generations, int steps_per_generation,
            int delay, double fraction, int seed, double temp)
{
    // Get variable screen information
    struct fb_var_screeninfo vinfo;
    if (ioctl(fbfd, FBIOGET_VSCREENINFO, &vinfo) == -1) {
        exit(3);
    }
    
    if (vinfo.bits_per_pixel != 32) {
        exit(5);
    }

    World m(vinfo.yres, vinfo.xres, temp, seed);
    m.init(fraction);

    long screensize = vinfo.xres * vinfo.yres * 4;

    // Map the device to memory
    uint32_t* fbp = (uint32_t*)mmap(0, screensize, PROT_READ | PROT_WRITE, MAP_SHARED, fbfd, 0);
    if (long(fbp) == -1) {
        exit(4);
    }

    uint32_t green = 0x0000ff00; // aarrggbb, assuming bits_per_pixel == 32
    uint32_t red = 0x00ff0000; 

    // write matrix to framebuffer
    auto setter = [&green, &red](auto x){return x == 1? green : red;};
    transform(begin(m.data()), end(m.data()), fbp, setter);
    
    for (int i = 0; i < generations; i++)
    {
        this_thread::sleep_for(chrono::milliseconds(delay));
        m.update(steps_per_generation);
        transform(begin(m.data()), end(m.data()), fbp, setter);
        msync(fbp, screensize, MS_SYNC);
    }
        
    // cleanup
    munmap(fbp, screensize);
    close(fbfd);
    return 0;
}


int main(int argc, char* argv[])
{
    if (argc > 1 and argv[1][0] == '?')
    {
        printf("Usage: "
               "%s <generations> <steps_per_generation> <delay (ms)> <init fraction> <seed> <prefer_txt> <temp>\n",
               argv[0]);
        exit(0);
    }
    
    int generations = (argc > 1) ? atoi(argv[1]) : 10000;
    int steps_per_generation = (argc > 2) ? atoi(argv[2]) : 1000;
    int delay = (argc > 3) ? atoi(argv[3]) : 200;
    double fraction = (argc > 4) ? atof(argv[4]) : 0.3;
    int seed = (argc > 5) ? atoi(argv[5]) : 0;
    bool prefer_txt = (argc > 6) ? bool(atoi(argv[6])) : false;
    double temp = (argc > 7) ? atof(argv[7]) : 1.0;
    
    if (not prefer_txt)
    {
        // Try to open the framebuffer for reading and writing
        int fbfd = open("/dev/fb0", O_RDWR);
        if (fbfd != -1) {
            return main_fb(fbfd, generations, steps_per_generation, delay, fraction, seed, temp);
        }
    }
    return main_txt(generations, steps_per_generation, delay, fraction, seed, temp);
}

