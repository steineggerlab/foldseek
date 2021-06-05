// Copyright 2019 Global Phasing Ltd.
//
// Span - span of array or std::vector.
// MutableVectorSpan - span of std::vector with insert() and erase()

#ifndef GEMMI_SPAN_HPP_
#define GEMMI_SPAN_HPP_

#include <algorithm>    // for find_if, find_if_not
#include <vector>
#include <stdexcept>    // for out_of_range
#include <type_traits>  // for remove_cv, conditional, is_const

namespace gemmi {

template<typename Item> struct MutableVectorSpan;

// Minimalistic Span, somewhat similar to C++20 std::span.
template<typename Item> struct Span {
  using iterator = Item*;
  using const_iterator = Item const*;
  using element_type = Item;
  using value_type = typename std::remove_cv<Item>::type;

  friend Span<const value_type>;
  friend MutableVectorSpan<value_type>;

  Span() = default;
  Span(iterator begin, std::size_t n) : begin_(begin), size_(n) {}

#if !defined(_MSC_VER) || _MSC_VER-0 >= 1926
  // constructor only for const Item, to allow non-const -> const conversion
  template<typename T=Item>
  Span(const Span<value_type>& o,
       typename std::enable_if<std::is_const<T>::value>::type* = 0)
#else
  // older MSVC don't like the version above
  Span(const Span<value_type>& o)
#endif
    : begin_(o.begin_), size_(o.size_) {}

  void set_begin(iterator begin) { begin_ = begin; }
  void set_size(std::size_t n) { size_ = n; }
  const_iterator begin() const { return begin_; }
  const_iterator end() const { return begin_ + size_; }
  iterator begin() { return begin_; }
  iterator end() { return begin_ + size_; }

  Item& front() { return *begin_; }
  const Item& front() const { return *begin_; }
  Item& back() { return *(begin_ + size_ - 1); }
  const Item& back() const { return *(begin_ + size_ - 1); }

  const Item& operator[](std::size_t i) const { return *(begin_ + i); }
  Item& operator[](std::size_t i) { return *(begin_ + i); }

  Item& at(std::size_t i) {
    if (i >= size())
      throw std::out_of_range("item index ouf of range: #" + std::to_string(i));
    return *(begin_ + i);
  }
  const Item& at(std::size_t i) const {
    return const_cast<Span*>(this)->at(i);
  }

  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  explicit operator bool() const { return size_ != 0; }

  template<typename Iter> Span<Item> sub(Iter first, Iter last) {
    return Span<Item>(&*first, last - first);
  }

  template<typename F, typename V=Item> Span<V> subspan(F&& func) {
    auto group_begin = std::find_if(this->begin(), this->end(), func);
    auto group_end = std::find_if_not(group_begin, this->end(), func);
    return Span<V>(&*group_begin, group_end - group_begin);
  }
  template<typename F> Span<const value_type> subspan(F&& func) const {
    using V = const value_type;
    return const_cast<Span*>(this)->subspan<F, V>(std::forward<F>(func));
  }

private:
  iterator begin_;
  std::size_t size_ = 0;
};

// Span of std::vector, implements insert() and erase().
template<typename Item> struct MutableVectorSpan : Span<Item> {
  using vector_type = std::vector<typename Span<Item>::value_type>;
  using iterator = typename Span<Item>::iterator;
  //friend Span<const value_type>;
  MutableVectorSpan() = default;
  MutableVectorSpan(Span<Item>&& p, vector_type* v)
    : Span<Item>(p), vector_(v) {}
  MutableVectorSpan(vector_type& v, iterator begin, std::size_t n)
    : Span<Item>(begin, n), vector_(&v) {}

  template<typename Iter> MutableVectorSpan<Item> sub(Iter first, Iter last) {
    return {Span<Item>::sub(first, last), vector_};
  }

  template<typename F> MutableVectorSpan<Item> subspan(F&& func) {
    return {Span<Item>::subspan(std::forward<F>(func)), vector_};
  }
  template<typename F> MutableVectorSpan<const Item> subspan(F&& func) const {
    return {Span<const Item>::subspan(std::forward<F>(func)), vector_};
  }

  iterator insert(iterator pos, Item&& item) {
    auto offset = this->begin_ - this->vector_->data();
    auto iter = vector_->begin() + (pos - this->vector_->data());
    auto ret = vector_->insert(iter, std::move(item));
    this->begin_ = vector_->data() + offset;
    ++this->size_;
    return &*ret;
  }

  void erase(iterator pos) {
    vector_->erase(vector_->begin() + (pos - vector_->data()));
    --this->size_;
  }

  bool is_beginning() const { return this->begin() == vector_->data(); }
  bool is_ending() const { return this->end() == vector_->data() + vector_->size(); }

private:
  vector_type* vector_ = nullptr;  // for insert() and erase()
};

} // namespace gemmi
#endif
