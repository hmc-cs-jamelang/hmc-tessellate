// Example program
#include <iostream>
#include <string>
#include <vector>
#include <iterator>

template <typename Iterator>
struct Sequence {
    using iterator = Iterator;
    using const_iterator = Iterator;

    Iterator begin_, end_;
    Sequence(Iterator begin, Iterator end)
        : begin_(begin), end_(end)
    { /* Done */ }

    Iterator begin() { return begin_; }
    Iterator begin() const { return begin_; }

    Iterator end() { return end_; }
    Iterator end() const { return end_; }
};

template <typename Iterator>
Sequence<Iterator> sequence(Iterator begin, Iterator end)
{
    return Sequence<Iterator> {begin, end};
}

template <typename BaseIterator, typename Predicate>
struct FilterIterator : BaseIterator {
private:
    const BaseIterator end_;
    const Predicate& predicate_;

    FilterIterator(BaseIterator start, BaseIterator end, Predicate predicate)
        : BaseIterator(start), end_(end), predicate_(predicate)
    {
        while (*this != end_ && !predicate_(**this)) {
            BaseIterator::operator++();
        }
    }

    FilterIterator(BaseIterator end, Predicate predicate)
        : BaseIterator(end), end_(end), predicate_(predicate)
    {}

public:
    template <typename Collection>
    static Sequence<FilterIterator> filter(const Collection& c, Predicate p)
    {
        return sequence(FilterIterator {begin(c), end(c), p},
                        FilterIterator {end(c), p});
    }

    FilterIterator& operator++()
    {
        do {
            BaseIterator::operator++();
        } while (*this != end_ && !predicate_(**this));
        return *this;
    }
};


using std::begin;
using std::end;

using FilterIterator::filter;

int main()
{
  std::vector<int> x;
  x.push_back(1);
  x.push_back(2);
  x.push_back(3);

//   for (auto y : x) { std::cout << y << std::endl; }
//   for (auto y : sequence(x.begin(), x.end())) { std::cout << y << std::endl; }
//   for (auto y : sequence(x.rbegin(), x.rend())) { std::cout << y << std::endl; }



  for (auto y : filter(x, [](int x) {return x > 1;})) {
      std::cout << y << std::endl;
  }
  return 0;
}
