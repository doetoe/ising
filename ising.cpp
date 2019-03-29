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
              generator_(), generator2_(),
              row_picker_(0, rows - 1), col_picker_(0, cols - 1), dist_()
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

    double get_temp() const
    {
        return 1./beta_;
    }
    
    // Average magnetization
    double net_magnetization() const
    {
        // cannot use data().sum(), because the data type cannot hold the sum in general
        return accumulate(begin(data()), end(data()), 0.) / // default operator is plus<T>
            (getRows() * getCols());
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
    
    void wolff_update(int n=1) 
    {
        // to be implemented
    }
    
    void print(const string& info) const
    {
        uint32_t R = getRows();
        uint32_t C = getCols();
        uint32_t startrow = 0;
        if (info != "")
        {
            startrow = 1;
            cout << info << endl;
        }
        
        for (uint32_t r = startrow; r < R; r++)
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
#include <sstream>
#include <vector>
#include <sys/ioctl.h>
#include <unistd.h>
#include <fcntl.h>
#include <linux/fb.h>
#include <sys/mman.h>

class Interaction
{
    World* world_;
    unsigned char c;
    bool show_info_;
    double delay_;
    double steps_per_generation_;
    
public:
    enum KeyAction {EXIT, CONTINUE};
        
    Interaction(World* world, double delay, uint32_t steps_per_generation)
            : world_(world), c('\0'), show_info_(false), delay_(delay),
              steps_per_generation_(steps_per_generation)
    {
        fcntl(STDIN_FILENO, F_SETFL, O_NONBLOCK);  // make the reads non-blocking
    }

    void change_temp(double factor)
    {
        world_->set_temp(world_->get_temp() * factor);
    }

    void change_delay(double factor)
    {
        delay_ *= factor;
    }

    void change_steps_per_generation(double factor)
    {
        steps_per_generation_ *= factor;
    }

    // Delay in integer milliseconds
    uint32_t get_delay() const
    {
        return uint32_t(delay_ + .5);
    }

    uint32_t get_steps_per_generation() const
    {
        return uint32_t(steps_per_generation_ + 0.5);
    }

    string info_string() const // could add help
    {
        //ostringstream ret;
        
        if (show_info_)
        {
            /*
            ret << "Temperature: " << world_->get_temp() << "  "
                << "Magnetization: " << world_->net_magnetization() << "  "
                << "Delay: " << get_delay() << " ms  "
                << "Steps per generation: " << get_steps_per_generation() << "  "
                << "Commands: hcfsmliwq";
            */
            auto format =
                "  Temperature: %.6f"
                "  Magnetization: %.3f"
                "  Delay: %d ms"
                "  Steps per generation: %d"
                "  Commands: hcfsmliwq";
            int len = snprintf(nullptr, 0, format,
                               world_->get_temp(), world_->net_magnetization(),
                               get_delay(), get_steps_per_generation()) + 1;
            vector<char> chars(len);
            snprintf(&chars[0], chars.size(), format,
                     world_->get_temp(), world_->net_magnetization(),
                     get_delay(), get_steps_per_generation());
            
            return string(&chars[0]);
        }
        else
        {
            return "";
        }
        // return ret.str();
    }

    void toggle_info()
    {
        show_info_ = !show_info_;
    }
    
    // Returns whether to exit
    KeyAction check_for_key()
    {
        if (read(fileno(stdin), &c, 1) > 0)
            // c = getchar();
        {
            switch (c)
            {
                case 'h': // hotter
                    change_temp(1.1);
                    break;
                case 'c': // colder
                    change_temp(1/1.1);
                    break;
                case 'f': // faster
                    change_delay(1/1.1);
                    break;
                case 's': // slower
                    change_delay(1.1);
                    break;
                case 'm': // more
                    change_steps_per_generation(1.1);
                    break;
                case 'l': // less
                    change_steps_per_generation(1/1.1);
                    break;
                case 'i': // info
                    toggle_info();
                    break;
                case 'w': // Wolff
                    world_->wolff_update();
                    break;
                case 'q':
                    return EXIT;
            }
        }
        return CONTINUE;
    }
};


int main_txt(int steps_per_generation,
             double delay, double fraction, int seed, double temp)
{
    struct winsize size;
    ioctl(STDOUT_FILENO,TIOCGWINSZ,&size);

    World m(size.ws_row, size.ws_col, temp, seed);
    // World m(5,13, temp);
    // printf("Size: %d x %d\n", size.ws_row, size.ws_col);
    // exit(0);
    m.init(fraction);
    m.print("");

    Interaction interaction(&m, delay, steps_per_generation);
    
    while (true)
    {
        if (interaction.check_for_key() == Interaction::EXIT)
        {
            break;
        }
        
        this_thread::sleep_for(std::chrono::milliseconds(interaction.get_delay()));
        m.update(interaction.get_steps_per_generation());
        m.print(interaction.info_string());
    }
    return 0;
}

// fbfd is the open (R/W) file descriptor of the framebuffer
int main_fb(int fbfd, int steps_per_generation,
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
    
    while (true)
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
    if ((argc > 1 and argv[1][0] == 'h') or argc > 7)
    {
        printf("Usage: "
               "%s <temp> <steps_per_generation> <delay (ms)> <init fraction> <seed> <prefer_txt>\n",
               argv[0]);
        exit(0);
    }
    
    double temp = (argc > 1) ? atof(argv[1]) : 1.0;
    int steps_per_generation = (argc > 2) ? atoi(argv[2]) : 1000;
    int delay = (argc > 3) ? atoi(argv[3]) : 200;
    double fraction = (argc > 4) ? atof(argv[4]) : 0.5;
    int seed = (argc > 5) ? atoi(argv[5]) : 0;
    bool prefer_txt = (argc > 6) ? bool(atoi(argv[6])) : false;
    
    if (not prefer_txt)
    {
        // Try to open the framebuffer for reading and writing
        int fbfd = open("/dev/fb0", O_RDWR);
        if (fbfd != -1) {
            return main_fb(fbfd, steps_per_generation, delay, fraction, seed, temp);
        }
    }
    return main_txt(steps_per_generation, delay, fraction, seed, temp);
}

