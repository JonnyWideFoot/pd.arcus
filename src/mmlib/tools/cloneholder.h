#ifndef __CLONEHOLDER
#define __CLONEHOLDER

//-------------------------------------------------
//
/// \brief  CloneHolder - holds classes which support clone(), simplifying their management.
///
/// \details 
/// A one member holder that can hold an object which supports clone(). 
/// An instance of this class should be owned by a parent class that requires a memory
/// managed pointer, but doesnt want to implement the full set of destructors and copy 
/// constructors to manage the 'new' allocated memory. This class then deals with any  
/// memory management issues, meaning that the owning class doesnt have to.
///
/// static_cast<T*> oddly seems to be required when not holding the direct base class, even though
/// the compiler should be able to assert this from T, oddness...
///
/// \author  Jon Rea 
template<typename T>
class CloneHolder
{
public:
	CloneHolder() : m_Data(NULL) {}
	~CloneHolder() { delete m_Data; }
	CloneHolder( const CloneHolder<T>& _clone ) { m_Data = static_cast<T*>(_clone.m_Data->clone()); }
	CloneHolder( const T& _clone ) { m_Data = static_cast<T*>(_clone.clone()); }
	CloneHolder( const T* _clone ) { m_Data = static_cast<T*>(_clone->clone()); }
	CloneHolder& operator=( const CloneHolder<T>& _clone ) { delete m_Data; m_Data = static_cast<T*>(_clone.m_Data->clone()); return *this; }
	CloneHolder& operator=( const T& _clone ) { delete m_Data; m_Data = static_cast<T*>(_clone.clone()); return *this; }
	CloneHolder& operator=( const T* _clone ) { delete m_Data; m_Data = static_cast<T*>(_clone->clone()); return *this; }
	T* operator->() { return m_Data; }
	const T* operator->() const { return m_Data; }
	const T* getPtr() const { return m_Data; }
	const T& data() const { return *m_Data; }
	T& data() { return *m_Data; }
	void setTo( const T& _clone ) { delete m_Data; m_Data = static_cast<T*>(_clone.clone()); }
	bool assigned() const { return m_Data != NULL; }
	void release() { delete m_Data; m_Data = NULL; return; }
	
private:
	T* m_Data;
};

#endif

