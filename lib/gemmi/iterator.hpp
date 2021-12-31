// Copyright 2018 Global Phasing Ltd.
//
// Iterators. Currently each of them is a BidirectionalIterator.

#ifndef GEMMI_ITERATOR_HPP_
#define GEMMI_ITERATOR_HPP_
#include <iterator>     // for bidirectional_iterator_tag
#include <type_traits>  // for remove_cv
#include <vector>

#ifdef  __INTEL_COMPILER
// warning #597: "X<T>::operator X<T>() const" will not be called for implicit
// or explicit conversions. That warning is triggered when templates
// StrideIter, IndirectIter and others are expanded with const Value.
# pragma warning disable 597
#endif

namespace gemmi {

// implements concept BidirectionalIterator
template <typename Policy>
struct BidirIterator : Policy {
  using value_type = typename std::remove_cv<typename Policy::value_type>::type;
  using difference_type = std::ptrdiff_t;
  using pointer = typename Policy::value_type*;
  using reference = typename Policy::reference;
  using iterator_category = std::bidirectional_iterator_tag;

  BidirIterator() = default;
  BidirIterator(Policy&& p) : Policy(p) {}

  BidirIterator& operator++() { Policy::increment(); return *this; }
  BidirIterator operator++(int) { BidirIterator x = *this; ++*this; return x; }
  BidirIterator& operator--() { Policy::decrement(); return *this; }
  BidirIterator operator--(int) { BidirIterator x = *this; --*this; return x; }
  bool operator==(const BidirIterator &o) const { return Policy::equal(o); }
  bool operator!=(const BidirIterator &o) const { return !Policy::equal(o); }
  reference operator*() { return Policy::dereference(); }
  pointer operator->() { return &Policy::dereference(); }
  using const_variant = BidirIterator<typename Policy::const_policy>;
  operator const_variant() const {
    return const_variant(static_cast<const Policy&>(*this));
  }
};

template<typename Value>
class StrideIterPolicy {
public:
  using value_type = Value;
  using reference = Value&;
  StrideIterPolicy() : cur_(nullptr), offset_(0), stride_(0) {}
  StrideIterPolicy(Value* ptr, std::size_t offset, size_t stride)
    : cur_(ptr), offset_(offset), stride_((unsigned)stride) {}
  void increment() { cur_ += stride_; }
  void decrement() { cur_ -= stride_; }
  bool equal(const StrideIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return cur_[offset_]; }
  using const_policy = StrideIterPolicy<Value const>;
  operator const_policy() const { return const_policy(cur_, offset_, stride_); }
private:
  Value* cur_;
  std::size_t offset_;
  unsigned stride_;
};
template<typename Value>
using StrideIter = BidirIterator<StrideIterPolicy<Value>>;


template<typename Redirect, typename Value>
class IndirectIterPolicy {
public:
  using value_type = Value;
  using reference = Value&;
  IndirectIterPolicy() : redir_(nullptr) {}
  IndirectIterPolicy(Redirect* redir, std::vector<int>::const_iterator cur)
    : redir_(redir), cur_(cur) {}
  void increment() { ++cur_; }
  void decrement() { --cur_; }
  bool equal(const IndirectIterPolicy& o) const { return cur_ == o.cur_; }
  Value& dereference() { return redir_->value_at(*cur_); }
  using const_policy = IndirectIterPolicy<Redirect const, Value const>;
  operator const_policy() const { return const_policy(redir_, cur_); }
  // TODO: what should be done with absent optional tags (*cur_ < 0)?
private:
  Redirect* redir_;
  std::vector<int>::const_iterator cur_; // points into positions
};
template<typename Redirect, typename Value>
using IndirectIter = BidirIterator<IndirectIterPolicy<Redirect, Value>>;


template<typename Vector, typename Value>
class UniqIterPolicy {
public:
  using value_type = Value;
  using reference = Value&;
  UniqIterPolicy() : vec_(nullptr), pos_(0) {}
  UniqIterPolicy(Vector* vec, std::size_t pos) : vec_(vec), pos_(pos) {}
  void increment() {
    // move to the first element of the next group
    const auto& key = (*vec_)[pos_].group_key();
    ++pos_;
    while (pos_ != vec_->size() && (*vec_)[pos_].group_key() == key)
      ++pos_;
  }
  void decrement() {
    --pos_; // now we are at the last element of the previous group
    const auto& key = (*vec_)[pos_].group_key();
    while (pos_ != 0 && (*vec_)[pos_-1].group_key() == key)
      --pos_; // move to the group beginning
  }
  bool equal(const UniqIterPolicy& o) const { return pos_ == o.pos_; }
  Value& dereference() { return (*vec_)[pos_]; }
  using const_policy = UniqIterPolicy<Vector const, Value const>;
  operator const_policy() const { return const_policy(vec_, pos_); }
private:
  Vector* vec_;
  std::size_t pos_;
};
template<typename Vector, typename Value>
using UniqIter = BidirIterator<UniqIterPolicy<Vector, Value>>;

template<typename Value, typename Vector=std::vector<Value>>
struct UniqProxy {
  Vector& vec;
  using iterator = UniqIter<Vector, Value>;
  iterator begin() { return {{&vec, 0}}; }
  iterator end() { return {{&vec, vec.size()}}; }
};
template<typename Value, typename Vector=std::vector<Value>>
struct ConstUniqProxy {
  const Vector& vec;
  using iterator = UniqIter<const Vector, const Value>;
  iterator begin() const { return {{&vec, 0}}; }
  iterator end() const { return {{&vec, vec.size()}}; }
};


template<typename Vector, typename Value>
class GroupingIterPolicy {
public:
  using value_type = Value;
  using reference = Value&;
  GroupingIterPolicy() = default;
  GroupingIterPolicy(const Value& span) : span_(span) {}
  void increment() {
    span_.set_begin(span_.end());
    span_.set_size(0);
    while (!span_.is_ending() &&
           span_.begin()->group_key() == span_.end()->group_key())
      span_.set_size(span_.size() + 1);
  }
  void decrement() {
    span_.set_begin(span_.begin() - 1);
    span_.set_size(1);
    while (!span_.is_beginning() &&
           span_.begin()->group_key() == (span_.begin() - 1)->group_key()) {
      span_.set_begin(span_.begin() - 1);
      span_.set_size(span_.size() + 1);
    }
  }
  bool equal(const GroupingIterPolicy& o) const {
    return span_.begin() == o.span_.begin();
  }
  Value& dereference() { return span_; }
  using const_policy = GroupingIterPolicy<Vector const, Value const>;
  operator const_policy() const { return const_policy(span_); }
private:
  Value span_;
};
template<typename Vector, typename Value>
using GroupingIter = BidirIterator<GroupingIterPolicy<Vector, Value>>;


template<typename Filter, typename Vector, typename Value>
class FilterIterPolicy {
public:
  using value_type = Value;
  using reference = Value&;
  FilterIterPolicy() : vec_(nullptr), pos_(0) {}
  FilterIterPolicy(const Filter* filter, Vector* vec, std::size_t pos)
      : filter_(filter), vec_(vec), pos_(pos) {
    while (pos_ != vec_->size() && !matches(pos_))
      ++pos_;
  }
  bool matches(std::size_t p) const { return filter_->matches((*vec_)[p]); }
  void increment() { while (++pos_ < vec_->size() && !matches(pos_)) {} }
  void decrement() { while (pos_ != 0 && !matches(--pos_)) {} }
  bool equal(const FilterIterPolicy& o) const { return pos_ == o.pos_; }
  Value& dereference() { return (*vec_)[pos_]; }
  using const_policy = FilterIterPolicy<Filter, Vector const, Value const>;
  operator const_policy() const { return const_policy(vec_, pos_); }
private:
  const Filter* filter_;
  Vector* vec_;
  std::size_t pos_;
};
template<typename Filter, typename Vector, typename Value>
using FilterIter = BidirIterator<FilterIterPolicy<Filter, Vector, Value>>;

template<typename Filter, typename Value>
struct FilterProxy {
  const Filter& filter;
  std::vector<Value>& vec;
  using iterator = FilterIter<Filter, std::vector<Value>, Value>;
  iterator begin() { return {{&filter, &vec, 0}}; }
  iterator end() { return {{&filter, &vec, vec.size()}}; }
};

template<typename Filter, typename Value>
struct ConstFilterProxy {
  const Filter& filter;
  const std::vector<Value>& vec;
  using iterator = FilterIter<Filter, const std::vector<Value>, const Value>;
  iterator begin() const { return {{&filter, &vec, 0}}; }
  iterator end() const { return {{&filter, &vec, vec.size()}}; }
};


template<typename Item>
struct ItemGroup {
  using element_type = Item;

  ItemGroup(Item* start, const Item* end)
      : size_(int(end - start)), extent_(int(end - start)), start_(start) {
    for (const Item* i = start + 1; i != end; ++i)
      if (i->group_key() != start->group_key())
        --size_;
  }

  struct iterator {
    Item* ptr;
    const Item* end;
    bool operator==(const iterator& o) const { return ptr == o.ptr; }
    bool operator!=(const iterator& o) const { return ptr != o.ptr; }
    iterator& operator++() {
      const Item* prev = ptr++;
      while (ptr != end && ptr->group_key() != prev->group_key())
        ++ptr;
      return *this;
    }
    Item& operator*() { return *ptr; }
    Item* operator->() { return ptr; }
  };
  iterator begin() { return iterator{start_, start_+extent_}; }
  iterator end() { return iterator{start_+extent_, start_+extent_}; }

  size_t size() const { return (size_t) size_; }
  int extent() const { return extent_; }
  bool empty() const { return size_ == 0; }
  Item& front() { return *start_; }
  const Item& front() const { return *start_; }
  Item& back() { return start_[extent_ - 1]; }
  const Item& back() const { return start_[extent_ - 1]; }

  // constant time unless sparse (extend_ > size_)
  Item& operator[](std::size_t i) {
    if (size_ == extent_ || i == 0)
      return start_[i];
    for (Item* ptr = start_ + 1; ; ++ptr)
      if (ptr->group_key() == start_->group_key())
        if (--i == 0)
          return *ptr;
  }
  const Item& operator[](std::size_t i) const {
    return const_cast<ItemGroup*>(this)->operator[](i);
  }

private:
  int size_ = 0;
  int extent_ = 0;
  Item* start_ = nullptr;
};

} // namespace gemmi
#endif
