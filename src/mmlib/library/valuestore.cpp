#include "global.h"

#include "valuestore.h"

namespace Library{

ValueStore* ValueStore::getSingleton()
{
	static ValueStore inst;
	return &inst;
}

ValueStore::ValueStore()
{
}
	


size_t ValueStore::addNewData(
	const std::string &_newname, 
	const std::string &_newdata
){
	int index = findName(_newname);
	if(index<0){
		// add name if it's not already there
		m_Name.push_back(_newname);
		index = m_Name.size() - 1;
	}

	ValueData newdata;
	newdata.m_NameIndex = index;
	newdata.m_Data      = _newdata;	
	m_Value.push_back(newdata);
	return m_Value.size() - 1;
}

bool ValueStore::hasName(size_t id, const std::string &_qname){
	if( id >= m_Value.size() ){
		throw(ArgumentException("Data ID given to ValueStore is invalid. No such piece of data."));
	}
	if( cmpstring( m_Name[ m_Value[id].m_NameIndex ] , _qname ) ) return true;
	return false;
}

const std::string &ValueStore::getRawData(size_t id) const
{
	if( id >= m_Value.size() ){
		throw(ArgumentException("Data ID given to ValueStore is invalid. No such piece of data."));
	}
	return m_Value[id].m_Data ;
}

void ValueStore::setRawData(size_t id, const std::string& _newData) 
{
	if( id >= m_Value.size() ){
		throw(ArgumentException("Data ID given to ValueStore is invalid. No such piece of data."));
	}
	m_Value[id].m_Data = _newData;
}

int ValueStore::findName(const std::string &_qname)
{
	for(size_t i=0; i < m_Name.size(); i ++ ){
		if(cmpstring( _qname, m_Name[i] ) ) return i;
	}
	return -1;
}



}
