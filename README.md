# The Ising Model #

This is a fast implementation of a simulation of spin dynamics according to the [Ising Model](https://en.wikipedia.org/wiki/Ising_model) on Linux. When available, e.g. when executing the program in a non-graphical mode (Ctrl-Alt-F1), the full framebuffer is used so that every pixel is a cell. When the framebuffer is not available, the present terminal window is used, in which every character is a cell. In both cases we have periodic boundaries.

Note that you may have to be a superuser to have permissions to use the framebuffer.

Finally, the text-based output will look better or worse depending on the terminal emulator you are using. The fastest and best looking one that I tried was `rxvt-unicode` (`urxvt`).

### Usage ###

Compile the executable by invoking `make`. This will generate the program `ising`. When run with a single argument `?`, usage information is displayed, namely 

    Usage: ./ising <generations> <steps_per_generation> <delay (ms)> <init fraction> <seed> <prefer_txt> <temp>

* <generations> (default 100000) is the number of generations you want to visualize
* <steps_per_generation> (default 1000) is the number of changes tried per generation
* <delay> (default 200 ms) is the time in milliseconds that you want to wait between generations. You can make this as short as you want, but when it is too short to finish drawing one generation before refreshing it, it won't look very smooth.
* <init fraction> is a number between 0 and 1 with the proportion of live cells
* <seed> is an integer that seeds the random generator
* set <prefer_txt> to 1 if you want to have text output, even if the framebuffer is available.
* <temp> is the temperature in natural units (k = 1). Values around 1 are interesting. At 3 it stays pretty random, at 0.3 it only goes in one direction everywhere.

When executed, it will visualize the specified number of generations, and return to the command line (note that the framebuffer contents will not be erased before it is explicitly overwritten). Ctrl-C to exit prematurely.

If you want timing information, compile `isingtime` using `make isingtime`. This can be executed like the other, but also as

    ./isingtime <generations> <rows> <cols>

in which case it will not display anything, only compute the specified number of generations on a world of the specified size, and output how long each update took in microseconds.

### Notes ###

* Most of the code was developed for the implementation of the [Game of Life](https://bitbucket.org/doetoe/life) automaton. 
* The execution in the framebuffer is visually very interesting
* Some window sizes don't work so well (only updates along the diagonal). Is this a bug?
* Right now it makes a number of updates for each generation, and redraws afterward. It might be interesting to only update the changed values.

### Contact ###

doetoe@protonmail.com
