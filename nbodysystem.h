#include <stdio.h>
#include <stdlib.h>

template <typename T>
class NBodySystem {
    public:
        NBodySystem(int nBodies) {}
        virtual ~NBodySystem() {}

        virtual void update(T delta) = 0;

    protected:
        NBodySystem() {}    // default constructor
        
        virtual void _init(int nBodies) = 0;
        virtual void _fin() = 0;
};
