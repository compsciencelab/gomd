# gomd

This is a test to write a scientific compute-intensive code using golang. Normally these codes are written in C/C++ but it would be nice to have some more modern alternative for rapid coding, yet mantaining the capabilities to decide data types for performance.

So far on the positive side

1. Nice, clean, short syntax C-like 
1. It is a typed language but it hides it reasonably well
1. Good compromise between readibility and coding restrictions
1. I like the auto indenting
1. I like the simple executable output for portability and distribution

On the negative side

1. There is nothing like numpy, scipy
1. There will never be until the language allows for operator overloading
1. It is as bad as not having even Sqrt for floats and two cast operations have to be used on top of the normal 64 bit version

In my opinion it could become the defacto standard for compute-intensive scientific computations if just function overloading and generic functions were introduced, but until then it is severely limited.
