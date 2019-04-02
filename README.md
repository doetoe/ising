# The Ising Model #

This is a fast implementation of a simulation of spin dynamics according to the [Ising Model](https://en.wikipedia.org/wiki/Ising_model) on Linux.

When available, e.g. when executing the program in a non-graphical mode (Ctrl-Alt-F1), the full framebuffer is used so that every pixel is a particle. When the framebuffer is not available, the present terminal window is used, in which every character is a particle. In both cases we have periodic boundaries.

Note that you may have to be a superuser to have permissions to use the framebuffer. To get a snapshot of the framebuffer, you can use a program like `fbgrab`.

Finally, the text-based output will look better or worse depending on the terminal emulator you are using. The fastest and best looking one that I tried was `rxvt-unicode` (`urxvt`).

### Usage ###

Compile the executable by invoking `make`. This will generate the program `ising`. When run with a single argument `?`, usage information is displayed, namely 

    Usage: ./ising <temp> <steps_per_generation> <delay (ms)> <init fraction> <seed> <prefer_txt>

* <temp> is the temperature in natural units (k = 1). The Curie temperature is around 2.269.
* <steps_per_generation> (default 1000) is the number of changes tried per generation
* <delay> (default 200 ms) is the time in milliseconds that you want to wait between generations. You can make this as short as you want, but when it is too short to finish drawing one generation before refreshing it, it won't look very smooth.
* <init fraction> is a number between 0 and 1 with the proportion of live cells
* <seed> is an integer that seeds the random generator
* set <prefer_txt> to 1 if you want to have text output, even if the framebuffer is available.

When executed, it will visualize the specified number of generations, and return to the command line (note that the framebuffer contents will not be erased before it is explicitly overwritten). Ctrl-C to exit prematurely.

In text mode, there is interaction as well. For now, you have to enter after your command(s). Commands are

* i     -- toggle info display
* h,c   -- hotter, colder
* f,s   -- faster, slower
* m,l   -- more, less
* w     -- step in Wolff cluster algorithm
* a     -- change algorithm (Wolff/Metropolis)
* q     -- quit

### Notes ###

* Most of the code was developed for the implementation of the [Game of Life](https://bitbucket.org/doetoe/life) automaton. 
* The execution in the framebuffer is visually very interesting
* Right now it makes a number of updates for each generation, and redraws afterward. It might be interesting to only update the changed values.

### The Ising Model ###

The Ising model is a simple model in statistical mechanics of a 2D lattice of magnatic dipoles, whose energy is obtained as the sum of the pairwise energies, where aligned direct (horizontal or vertical) neighbours contribute 1, and opposite dipoles contribute -1. It is known that there is a phase transition around a temperature of 2.269, the Curie temperature, above which it behaves as a paramagnetic material: there is no global magnetization, and below which it behaves as a ferromagnetic material, in which globally all spins tend to align.

In this implementation, successive states are generated following the Metropolis algorithm using the Boltzmann distribution for the specified temperature.

### Todo ###

* Make info display nices: flicker, layout, etc
* Update values not in display refresh loop
* Load/save states
* Roll (translate)
* Better help texts
* Other screen sizes
* Output/write statistics to file
* Display acceptance rate

### Contact ###

doetoe@protonmail.com
