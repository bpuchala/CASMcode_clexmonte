#ifndef CASM_clexmonte_misc_MapLike
#define CASM_clexmonte_misc_MapLike

namespace CASM {
namespace clexmonte {

struct MapLikeTraits {
  typedef std::size_t size_type;
  typedef std::size_t key_type;
};

template <typename T, bool IsConst>
class MapLikeIterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = std::conditional_t<IsConst, const T*, T*>;
  using reference = std::conditional_t<IsConst, const T&, T&>;

  using data_vector_ptr_type =
      std::conditional_t<IsConst, const std::vector<T>*, std::vector<T>*>;

  using key_type = MapLikeTraits::key_type;

  MapLikeIterator(data_vector_ptr_type data, const std::vector<bool>* assigned,
                  key_type key)
      : m_data(data), m_assigned(assigned), m_key(key) {
    skip_unassigned();
  }

  reference operator*() { return (*m_data)[m_key]; }
  pointer operator->() { return &(*m_data)[m_key]; }

  MapLikeIterator& operator++() {
    ++m_key;
    skip_unassigned();
    return *this;
  }

  MapLikeIterator operator++(int) {
    MapLikeIterator temp = *this;
    ++(*this);
    return temp;
  }

  bool operator==(const MapLikeIterator& other) const {
    return m_key == other.m_key;
  }
  bool operator!=(const MapLikeIterator& other) const {
    return !(*this == other);
  }

  key_type key() const { return m_key; }

 private:
  void skip_unassigned() {
    while (m_key < m_data->size() && !(*m_assigned)[m_key]) {
      ++m_key;
    }
  }

  data_vector_ptr_type m_data;
  const std::vector<bool>* m_assigned;
  key_type m_key;
};

/// \brief A data structure similar with an interface similar std::map but which
///     1) assigns a unique key instead of being provided one, and 2) can be
///     changed as determined necessary for efficiency of state-saving
///     operations.
///
/// \tparam T Must be default constructible
template <typename T>
class MapLike {
 public:
  using iterator = MapLikeIterator<T, false>;
  using const_iterator = MapLikeIterator<T, true>;
  typedef MapLikeTraits::size_type size_type;
  typedef MapLikeTraits::key_type key_type;

  MapLike() : m_n_assigned(0) {}

  void clear() {
    m_data.clear();
    m_assigned.clear();
    m_available.clear();
    m_n_assigned = 0;
  }

  /// \brief Get the number of assigned elements
  size_type size() const { return m_n_assigned; }

  /// \brief Get the number of allocated elements
  size_type size_allocated() const { return m_data.size(); }

  const_iterator begin() const {
    return const_iterator(&m_data, &m_assigned, 0);
  }

  const_iterator end() const {
    return const_iterator(&m_data, &m_assigned, m_data.size());
  }

  iterator begin() { return iterator(&m_data, &m_assigned, 0); }

  iterator end() { return iterator(&m_data, &m_assigned, m_data.size()); }

  const_iterator cbegin() const {
    return const_iterator(&m_data, &m_assigned, 0);
  }

  const_iterator cend() const {
    return const_iterator(&m_data, &m_assigned, m_data.size());
  }

  T& operator[](key_type key) { return m_data[key]; }

  T const& operator[](key_type key) const { return m_data[key]; }

  size_type count(key_type key) const {
    if (key < m_data.size()) {
      return m_assigned[key] ? 1 : 0;
    }
    return 0;
  }

  const_iterator find(key_type key) const {
    if (key < m_data.size() && m_assigned[key]) {
      return const_iterator(m_data, m_assigned, key);
    }
    return end();
  }

  /// \brief Insert an object and return the key assigned to it
  key_type insert() {
    key_type key;
    if (m_available.empty()) {
      key = m_data.size();
      m_data.emplace_back();
      m_assigned.push_back(true);
    } else {
      key = m_available.back();
      m_available.pop_back();
      m_assigned[key] = true;
    }
    m_n_assigned++;
    return key;
  }

  /// \brief Insert an object with a specific key if it is not already present
  ///
  /// \param key The key
  /// \return is_new_object, with the value `true` if a new object was inserted,
  ///     and the value `false `otherwise.
  bool insert(key_type key) {
    if (key < m_data.size()) {
      if (m_assigned[key]) {
        // Is allocated and already assigned - do nothing
        return false;

      } else {
        // Is allocated, but not already present - assign it

        // remove key from available:
        auto it = std::find(m_available.begin(), m_available.end(), key);
        if (it != m_available.end()) {
          // swap with the last element, then remove the last element
          std::iter_swap(it, m_available.end() - 1);
          m_available.pop_back();
        }
        m_assigned[key] = true;
        m_n_assigned++;
        return true;
      }
    } else {
      // Is not allocated, is not present - allocate it
      for (size_type i = m_data.size(); i < key; ++i) {
        m_available.push_back(i);
      }
      m_data.resize(key + 1);
      m_assigned.resize(key + 1, false);
      m_assigned[key] = true;
      m_n_assigned++;
      return true;
    }
  }

  /// \brief Erase an object with a specific key
  ///
  /// - This method will not erase the object from the vector, but will mark it
  ///   as unassigned and add it to the available keys.
  ///
  /// \param key The key of the object to erase
  void erase(key_type key) {
    m_n_assigned--;
    m_assigned[key] = false;
    m_available.push_back(key);
  }

 private:
  std::vector<T> m_data;

  std::vector<bool> m_assigned;

  std::vector<Index> m_available;

  size_type m_n_assigned;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
