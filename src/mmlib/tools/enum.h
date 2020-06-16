#ifndef __ENUM_TOOLS
#define __ENUM_TOOLS







//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details 
/// Code taken from:
/// http://www.codeproject.com/cpp/InheritEnum.asp
/// By Jon Rea of use in the MMLib project.
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
template <typename EnumT, typename BaseEnumT>
class InheritEnum
{
public:
  InheritEnum() {}
  InheritEnum(EnumT e)
	: enum_(e)
  {}

  InheritEnum(BaseEnumT e)
	: baseEnum_(e)
  {}

  explicit InheritEnum( int val )
	: enum_(static_cast<EnumT>(val))
  {}

  operator EnumT() const { return enum_; }
private:
  // Note - the value is declared as a union mainly for as a debugging aid. If
  // the union is undesired and you have other methods of debugging, change it
  // to either of EnumT and do a cast for the constructor that accepts BaseEnumT.
  union
  {
	EnumT enum_;
	BaseEnumT baseEnum_;
  };
};

#endif

