#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REFCOUNTEDPTR_H_
#define _REFCOUNTEDPTR_H_

#include <cstdlib>
#include <string>
#include <typeinfo>
#include "MayDay.H"
#include "Arena.H"
#include "CH_assert.H"
#ifdef DEBUGRCP
  #include "parstream.H"
#endif
#include "BaseNamespaceHeader.H"

#ifdef CH_USE_MEMORY_TRACKING
extern unsigned long long int ch_memcount;
#endif

#ifdef DEBUGRCP
  #define RCPDBG(x) x
#else
  #define RCPDBG(x) (void)0;
#endif

//--Traits

/*--------------------------------------------------------------------*
 * Provides a type with inverse constantness from T
 *--------------------------------------------------------------------*/

/// General case T is non-const and inverse is const
template <typename T>
struct RCPTypeTr
{
  typedef const T InvertConstType;
};

/// Specialization where T is const and inverse is non-const
template <typename T>
struct RCPTypeTr<const T>
{
  typedef T InvertConstType;
};

/*--------------------------------------------------------------------*
 * Checks if two types are the same
 *--------------------------------------------------------------------*/

/// General case T and S are not the same
template <typename T, typename S>
struct RCPPairTr
{
  enum
  {
    IsSame = 0
  };
};

/// Specialization for same type
template <typename T>
struct RCPPairTr<T, T>
{
  enum
  {
    IsSame = 1
  };
  typedef T SamePolicyType;
};

//--Policies

/*--------------------------------------------------------------------*
 * For destroying pointer types
 *--------------------------------------------------------------------*/

/// Policies for when T is a scalar object
struct RCPPointerPolicy
{
  template <typename T>
  static void destroy(T* a_pointer)
    {
      delete a_pointer;
    }
};

/// Policies for when T is an array
struct RCPArrayPolicy
{
  template <typename T>
  static void destroy(T* a_array)
    {
      delete[] a_array;
    }
};

/// Policies for when T is created via malloc (or otherwise requires free)
struct RCPFreePolicy
{
  template <typename T>
  static void destroy(T* a_mem)
    {
      free((void*)a_mem);
    }
};

/*--------------------------------------------------------------------*/
///  A reference-counting handle class.
/**

\tparam T   Type of object pointed to
\tparam OP  Object Policy for T.  If T* is an object pointer, use the
            RCPPointerPolicy.  If T* is an array, use the
            RCPArrayPolicy to get proper destruction.  If T* requires
            deletion with 'free', use the RCPFreePolicy.  By default
            RCPPointerPolicy is used.

This is to be used as a pointer to class T.  This will feel and smell
just like a built-in pointer except:

        -#  There is no need to call delete on the pointer.
        -#  The default copy constructor and assignment implement ref-counting.
        -#  The user may call isNonUnique to determine if this pointer is
            the only pointer to the data.  This can be used to hide the
            ref-counting behavior of a class.
        -#  Checks for dereference of a null pointer.

This class is completely inlined.

typical usage:

\code

{
  Box b;
  IntVect a;
  //refCount() == 1
  RefCountedPtr<IntVectSet> set(new IntVectSet());

  // still just one IntVectSet, but refCount()==2
  RefCountedPtr<IntVectSet>  otherSet(set);

  // Like a pointer, modifying one modifies the other
  otherSet->define(b);
  (*set)|=a;
  {
    RefCountedPtr<IntVectSet> anotherSet; // null
    anotherSet = otherSet ;  //now all three have refCount()==3
  }//anotherSet out of scope, so deleted.  IntVectSet NOT deleted.

  // set.refCount() == 2
  // otherSet.refCount() == 2;
  // otherset == set;

}
// when all RefCountedPtr's for a given object are deleted, the last
// one calls 'delete' on the member pointer.
\endcode

*//*------------------------------------------------------------------*/

template <typename T,
          typename OP = RCPPointerPolicy>  // OP means ObjectPolicy
class RefCountedPtr
{
public:
  // types
  typedef T value_type;
  typedef T* pointer;
  typedef T& reference;
  typedef RefCountedPtr<T, OP> Self;

  /// construct with POD pointer
  explicit RefCountedPtr(pointer const ptr = 0);

  /// copy
  RefCountedPtr(const Self& other);

  /// const converting copy adds const qualifier.
  /** the related converting copy could be used but this avoids the
      dynamic_cast when just adding qualifiers
  */
  RefCountedPtr(
    const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP>& other);

  /// related converting copy
  /** allows conversion between related POD pointers.  E.g. base to derived
      (requires dynamic_cast) or derived to base.
  */
  template <typename S>
  RefCountedPtr(const RefCountedPtr<S, OP>& other);

  // COMPILE TIME ERROR HANDLING
  // If we copy construct a RCP with a different object policy, the compiler
  // may try implicit conversion to POD pointer and then use the pointer
  // constructor.  This is bad since the pointer will then be delete twice.
  // These copy constructors are designed to capture that case (the compiler
  // prefers them over conversion) and force an error.
  template <typename OP2>
  RefCountedPtr(const RefCountedPtr<T, OP2>& other)
    {
      // A compiler error here means 'other' is not using the same object
      // policy as 'this'
      //typedef typename RCPPairTr<OP, OP2>::SamePolicyType Check;
      MayDay::Error("Invalid runtime code");
    }
  template <typename OP2>
  RefCountedPtr(
    const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP2>& other)
    {
      // A compiler error here means 'other' is not using the same object
      // policy as 'this'
      //typedef typename RCPPairTr<OP, OP2>::SamePolicyType Check;
      MayDay::Error("Invalid runtime code");
    }

  /// Destructor
  ~RefCountedPtr();

  /// assignment operator copies pointer member
  /** copies pointer member and integer pointer, decreases 'this' refcount
      before assignment, then increases refcount of this=rhs.
  */
  Self& operator=(const Self& rhs);

  /// const converting assignment operator adds const qualifier.
  Self& operator=(
    const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP>& rhs);

  /// related assignment operator
  /** allows conversion between related POD pointers.  E.g. base to derived
      (requires dynamic_cast) or derived to base.
  */
  template <typename S>
  Self& operator=(const RefCountedPtr<S, OP>& rhs);

  // COMPILE TIME ERROR HANDLING
  // Similar to copy except an error of ambiguous operator overloading is
  // seen instead of implicit conversion.  These routines make the error
  // more clear
  template <typename OP2>
  Self& operator=(const RefCountedPtr<T, OP2>& rhs)
    {
      // A compiler error here means 'rhs' is not using the same object
      // policy as 'this'
      //typedef typename RCPPairTr<OP, OP2>::SamePolicyType Check;
      MayDay::Error("Invalid runtime code");
    }
  template <typename OP2>
  Self& operator=(
    const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP2>& rhs)
    {
      // A compiler error here means 'rhs' is not using the same object
      // policy as 'this'
      //typedef typename RCPPairTr<OP, OP2>::SamePolicyType Check;
      MayDay::Error("Invalid runtime code");
    }

  /// dereference access operator used like a pointer derefence access function.
  pointer operator ->() const;

  /// Return a const reference to the pointer and don't ask questions
  /** Same as operator -> but obtains the reference even if the pointer has
   *  not yet been assigned.  Used in some early setup where we want to track
   *  the RCP but it hasn't been defined yet.  If you want the pointer
   *  and know it has been allocated, -> should be preferred.
   */
  pointer const& getRefToThePointer() const;

  /// pointer dereference.
  reference operator *() const;

  /// [] access
  reference operator [](const int a_idx) const;

  /// auto conversion to regular pointer where required.
  /** kind of dangerous.  Deleting the returned pointer will cause havoc.
   */
  operator T* () const;

  bool isNull() const;

  /// true when refcount is one.
  bool isNonUnique() const;

  int refCount() const;

  void swap(Self& b);

  /// Will never delete ptr_
  /**
     useful in functions that need to return a pointer which, under some
     circumstances, we want to own, and under others we don't.
  */
  void neverDelete();

  template<typename T2, typename OP2>
  friend class RefCountedPtr;

  static const std::string& name()
  {
    return name_;
  }

  // These do not increment and decrement the reference count.  Rather they
  // relate to memory tracking.  However, deletion of the pointer is in
  // decrement
  void increment();
  void decrement();
protected:
  pointer ptr_;
  int* refCount_;
  static std::string name_;
  static BArena s_Arena;
  static int size_;
private:

};

template<typename T, typename OP>
std::string RefCountedPtr<T, OP>::name_(
  std::string("RefCountedPtr ")+std::string(typeid(T).name()));

//template<typename T, typename OP>
//BArena RefCountedPtr<T>::s_Arena(name_.c_str());
template<typename T, typename OP>
BArena RefCountedPtr<T, OP>::s_Arena(typeid(T).name());

template<typename T, typename OP>
int RefCountedPtr<T, OP>::size_ = sizeof(T)+2*sizeof(int);

template <typename T, typename OP>
inline void
RefCountedPtr<T, OP>::increment()
{
#ifdef CH_USE_MEMORY_TRACKING
  ch_memcount+=size_;
  s_Arena.bytes += size_;
  refCount_[1] = size_;
  if (s_Arena.bytes > s_Arena.peak)
    {
      s_Arena.peak = s_Arena.bytes;
    }
#endif
}

template <typename T, typename OP>
inline
void RefCountedPtr<T, OP>::decrement()
{
#ifdef CH_USE_MEMORY_TRACKING
  ch_memcount-=refCount_[1];
  s_Arena.bytes -= refCount_[1];
#endif
  RCPDBG(pout() << "====> Destroying " << ptr_ << std::endl;)
  OP::destroy(ptr_);
  ptr_ = 0;
  delete [] refCount_;
  refCount_ = 0;
}

template <typename T, typename OP>
inline
RefCountedPtr<T, OP>::RefCountedPtr(pointer const ptr)
  : ptr_(ptr),
    refCount_(0)
{
  RCPDBG(pout() << "POD pointer constructor\n";)
  if (ptr_)
    {
      refCount_ = new int[2];
      if (refCount_ == 0)
        MayDay::Error("RefCountedPtr::RefCountedPtr(T* ptr) out of memory");
      *refCount_ = 1;
      increment();
      RCPDBG(pout() << "====> Creating " << ptr << std::endl;)
    }
}

template <typename T, typename OP>
inline
RefCountedPtr<T, OP>::RefCountedPtr(const Self& other)
  : ptr_(other.ptr_),
    refCount_(other.refCount_)
{
  RCPDBG(pout() << "standard copy " << ptr_ << std::endl;)
  if (refCount_ != 0)
    ++(*refCount_);
}

template <typename T, typename OP>
inline
RefCountedPtr<T, OP>::RefCountedPtr(
  const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP>& other)
  : ptr_(other.ptr_),
    refCount_(other.refCount_)
{
  RCPDBG(pout() << "const conversion copy " << ptr_ << std::endl;)
  if (refCount_ != 0)
    ++(*refCount_);
}

template <typename T, typename OP>
template <typename S>
inline
RefCountedPtr<T, OP>::RefCountedPtr(const RefCountedPtr<S, OP>& other)
  : ptr_(dynamic_cast<T*>(other.ptr_)),
    refCount_(other.refCount_)
{
  RCPDBG(pout() << "related conversion copy " << ptr_ << std::endl;)
  // Note: a failure here means you cannot convert S* into T* as they are
  //       not properly related.
  if (other.ptr_ != NULL) CH_assert(ptr_ != NULL);
  if (refCount_ != 0)
    ++(*refCount_);
}

template <typename T, typename OP>
inline
RefCountedPtr<T, OP>::~RefCountedPtr()
{
  if (refCount_ != 0 && --(*refCount_) == 0)
    {
      decrement();
    }
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::Self&
RefCountedPtr<T, OP>::operator =(const Self& rhs)
{
  RCPDBG(pout() << "standard assign " << rhs.ptr_ << std::endl;)
  if (ptr_ != rhs.ptr_)
    {
      if (refCount_ != 0 && --(*refCount_) == 0)
        {
          decrement();
        }
      ptr_ = rhs.ptr_;
      refCount_ = rhs.refCount_;
      if (refCount_ != 0)
        ++(*refCount_);
    }
  return *this;
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::Self&
RefCountedPtr<T, OP>::operator =(
  const RefCountedPtr<typename RCPTypeTr<T>::InvertConstType, OP>& rhs)
{
  RCPDBG(pout() << "const conversion assign " << rhs.ptr_ << std::endl;)
  if (ptr_ != rhs.ptr_)
    {
      if (refCount_ != 0 && --(*refCount_) == 0)
        {
          decrement();
        }
      ptr_ = rhs.ptr_;
      refCount_ = rhs.refCount_;
      if (refCount_ != 0)
        ++(*refCount_);
    }
  return *this;
}

template <typename T, typename OP>
template <typename S>
inline typename RefCountedPtr<T, OP>::Self&
RefCountedPtr<T, OP>::operator =(const RefCountedPtr<S, OP>& rhs)
{
  RCPDBG(pout() << "related conversion assign " << rhs.ptr_ << std::endl;)
  if (ptr_ != rhs.ptr_)
    {
      if (refCount_ != 0 && --(*refCount_) == 0)
        {
          decrement();
        }
      ptr_ = dynamic_cast<T*>(rhs.ptr_);
      // Note: a failure here means you cannot convert S* into T* as they are
      //       not properly related.
      if (rhs.ptr_ != NULL) CH_assert(ptr_ != NULL);
      refCount_ = rhs.refCount_;
      if (refCount_ != 0)
        ++(*refCount_);
    }
  return *this;
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::pointer
RefCountedPtr<T, OP>::operator ->() const
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator ->() on null pointer");
  return ptr_;
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::pointer const&
RefCountedPtr<T, OP>::getRefToThePointer() const
{
  return ptr_;
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::reference
RefCountedPtr<T, OP>::operator *() const
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator *() on null pointer");
  return *ptr_;
}

template <typename T, typename OP>
inline typename RefCountedPtr<T, OP>::reference
RefCountedPtr<T, OP>::operator [](const int a_idx) const
{
  CH_assert(ptr_);
  return ptr_[a_idx];
}

template <class T, typename OP>
inline
RefCountedPtr<T, OP>::operator T* () const
{
  RCPDBG(pout() << "auto conversion " << ptr_ << std::endl;)
  return ptr_;
}

template <typename T, typename OP>
inline bool
RefCountedPtr<T, OP>::isNull() const
{
  return (ptr_ == 0);
}

template <typename T, typename OP>
inline bool
RefCountedPtr<T, OP>::isNonUnique() const
{
  return refCount_ == 0 ? false : *refCount_ != 1;
}

template <typename T, typename OP>
inline int
RefCountedPtr<T, OP>::refCount() const
{
  return refCount_ == 0 ? 0 : *refCount_;
}

template <typename T, typename OP>
inline void
RefCountedPtr<T, OP>::swap(Self& b)
{
  pointer s = ptr_;
  ptr_ = b.ptr_;
  b.ptr_ = s;

  int* r = refCount_;
  refCount_ = b.refCount_;
  b.refCount_ = r;
}

template <typename T, typename OP>
inline void
RefCountedPtr<T, OP>::neverDelete()
{
  ++ *refCount_;
}

#include "BaseNamespaceFooter.H"
#endif
