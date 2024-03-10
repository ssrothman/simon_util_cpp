#ifndef SIMONTOOLS_RECURSIVE_REDUCE_H
#define SIMONTOOLS_RECURSIVE_REDUCE_H

template<typename T>
struct has_const_iterator
{
private:
    template<typename C> static char test(typename C::const_iterator*);
    template<typename C> static int  test(...);
public:
    enum { value = sizeof(test<T>(0)) == sizeof(char) };
};

template<class ValueType, class Function = std::plus<ValueType>>
auto recursive_reduce(const ValueType& input, ValueType init, const Function& f)
{
    return f(init, input);
}


template<class Container, class ValueType, class Function = std::plus<ValueType>, 
    typename = std::enable_if_t<has_const_iterator<Container>::value, void>>
auto recursive_reduce(const Container& input, ValueType init, const Function& f = std::plus<ValueType>())
{
    for (const auto &element: input) {
        auto result = recursive_reduce(element, ValueType{}, f);
        init = f(init, result);
    }

    return init;
}


#endif
