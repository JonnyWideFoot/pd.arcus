#ifndef __VALUESTORE_H
#define __VALUESTORE_H

namespace Library
{
	class ValueStore;






//-------------------------------------------------
//
/// \brief Class to store a given value in the ValueStore singleton 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API ValueData
	{
		friend class ValueStore;
	public:
		//setDataString(std::string _data){ m_Data = _data; }
		//const std::string &getDataString() const;
		//const float &getDataInteger() const;
		//const double &getDataDouble() const;
	private:
		int m_NameIndex;
		std::string m_Data;
	};






//-------------------------------------------------
//
/// \brief Central Value store which can be used to associate arbitrary data with names 
///
/// \details 
/// Singleton pattern - there's only one value store managing all the 
/// values. This avoids confusion about the indices handed out making them unique
/// The Value store acts a a giant warehouse for arbitrary strings which can be retrived
/// by their number whihc is created when a piece of data is added (like a catalogue
/// number.
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API ValueStore
	{
	public:
		static ValueStore* getSingleton(); /// 'Singleton Pattern' Implementation
		ValueStore();

		/// adds a new piece of data to the store returning the ID of the data
		/// which can be used to retrieve it.
		size_t addNewData(const std::string &_newname, const std::string &_newdata); 

		/// check if a certain piece of data has a certain name associated with it
		bool hasName(size_t id, const std::string &_qname);
		const std::string& getRawData(size_t id) const;
		void setRawData(size_t id, const std::string& _newData);

	private:
		/// Search for a name and return the index or -1 if the name is unique (new)
		int findName(const std::string &_qname);

		std::vector<ValueData> m_Value;  // store all the data (with name indices)
		std::vector<std::string> m_Name;  // store all the names
	};
}

#endif

