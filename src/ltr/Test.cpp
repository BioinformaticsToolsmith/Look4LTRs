#include <iostream>
using namespace std;
#include <memory>
 
class Rectangle {
    int length;
    int breadth;
 
public:
    Rectangle(int l, int b)
    {
        length = l;
        breadth = b;
    }

    ~Rectangle(){
        std::cout<< "Deleting ..." << std::endl;
    }
 
    int area()
    {
        return length * breadth;
    }
};
 
int main()
{
 
    auto ptr = new Rectangle(10, 5);
    shared_ptr<Rectangle> P1(ptr);
    // This'll print 50
    cout << P1->area() << endl;
 
    shared_ptr<Rectangle> P2;
    P2 = P1;
 
    // This'll print 50
    cout << P2->area() << endl;
 
    // This'll now not give an error,
    cout << P1->area() << endl;
 
    delete ptr;

    // This'll also print 50 now
    // This'll print 2 as Reference Counter is 2
    cout << P1.use_count() << endl;
    return 0;
}