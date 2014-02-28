/*
 * LMCParam.hh
 *
 *  Created on: Jan 14, 2014
 *      Author: nsoblath
 */

#ifndef LMCPARAM_HH_
#define LMCPARAM_HH_

#include "LMCException.hh"

#include "document.h"

#include <deque>
#include <map>
#include <sstream>
#include <string>

namespace locust
{
    class ParamValue;
    class ParamArray;
    class ParamNode;

    class Param
    {
        public:
            Param();
            Param(const Param& orig);
            virtual ~Param();

            virtual Param* Clone() const;

            virtual bool IsNull() const;
            virtual bool IsValue() const;
            virtual bool IsArray() const;
            virtual bool IsNode() const;

            ParamValue& AsValue();
            ParamArray& AsArray();
            ParamNode& AsNode();

            const ParamValue& AsValue() const;
            const ParamArray& AsArray() const;
            const ParamNode& AsNode() const;

            /// Assumes that the parameter is a value, and returns a reference to itself.
            const ParamValue& operator()() const;
            /// Assumes that the parameter is a value, and returns a reference to itself.
            ParamValue& operator()();

            /// Assumes that the parameter is an array, and returns a reference to the Param at aIndex.
            /// Behavior is undefined if aIndex is out-of-range.
            const Param& operator[]( unsigned aIndex ) const;
            /// Assumes that the parameter is an array, and returns a reference to the Param at aIndex.
            /// Behavior is undefined if aIndex is out-of-range.
            Param& operator[]( unsigned aIndex );

            /// Assumes that the parameter is a node, and returns a reference to the Param corresponding to aName.
            /// Throws an Exception if aName is not present.
            const Param& operator[]( const std::string& aName ) const;
            /// Assumes that the parameter is a node, and returns a reference to the Param corresponding to aName.
            /// Adds a new Value if aName is not present.
            Param& operator[]( const std::string& aName );

            virtual std::string ToString() const;

            static unsigned sIndentLevel;
    };

    class ParamValue : public Param
    {
        public:
            ParamValue();
            template< typename XStreamableType >
            ParamValue( XStreamableType aStreamable );
            ParamValue(const ParamValue& orig);
            virtual ~ParamValue();

            virtual Param* Clone() const;

            virtual bool IsNull() const;
            virtual bool IsValue() const;

            const std::string& Get() const;
            template< typename XValType >
            XValType Get() const;

            template< typename XStreamableType >
            ParamValue& operator<<( const XStreamableType& aStreamable );

            virtual std::string ToString() const;

        protected:
            std::string fValue;

    };

    template< typename XStreamableType >
    ParamValue::ParamValue( XStreamableType aStreamable ) :
            Param(),
            fValue()
    {
        (*this) << aStreamable;
    }

    template< typename XValType >
    XValType ParamValue::Get() const
    {
        XValType tReturn;
        std::stringstream tBuffer;
        tBuffer << fValue;
        tBuffer >> tReturn;
        return tReturn;
    }

    template< typename XStreamableType >
    ParamValue& ParamValue::operator<<( const XStreamableType& a_streamable )
    {
        std::stringstream tBuffer;
        tBuffer << a_streamable;
        fValue = tBuffer.str();
        return *this;
    }


    class ParamArray : public Param
    {
        public:
            typedef std::deque< Param* > contents;
            typedef contents::iterator iterator;
            typedef contents::const_iterator const_iterator;
            typedef contents::reverse_iterator reverse_iterator;
            typedef contents::const_reverse_iterator const_reverse_iterator;
            typedef contents::value_type contents_type;

        public:
            ParamArray();
            ParamArray( const ParamArray& orig );
            virtual ~ParamArray();

            virtual Param* Clone() const;

            virtual bool IsNull() const;
            virtual bool IsArray() const;

            unsigned Size() const;
            bool Empty() const;

            /// sets the size of the array
            /// if smaller than the current size, extra elements are deleted
            void Resize( unsigned a_size );

            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Throws a Exception if aName is not present or is not of type ParamValue
            const std::string& GetValue( unsigned aIndex ) const;
            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Throws a Exception if aName is not present or is not of type ParamValue
            template< typename XValType >
            XValType GetValue( unsigned aIndex ) const;

            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Returns aDefault if aName is not present or is not of type ParamValue
            const std::string& GetValue( unsigned aIndex, const std::string& aDefault ) const;
            const std::string& GetValue( unsigned aIndex, const char* aDefault ) const;
            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Returns aDefault if aName is not present or is not of type ParamValue
            template< typename XValType >
            XValType GetValue( unsigned aIndex, XValType aDefault ) const;

            /// Returns a pointer to the Param corresponding to aName.
            /// Returns NULL if aName is not present.
            const Param* At( unsigned aIndex ) const;
            /// Returns a pointer to the Param corresponding to aName.
            /// Returns NULL if aName is not present.
            Param* At( unsigned aIndex );

            /// Returns a pointer to the ParamValue (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            const ParamValue* ValueAt( unsigned aIndex ) const;
            /// Returns a pointer to the ParamValue (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            ParamValue* ValueAt( unsigned aIndex );
            /// Returns a pointer to the ParamValue (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            const ParamValue* ValueAtFast( unsigned aIndex ) const;
            /// Returns a pointer to the ParamValue (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            ParamValue* ValueAtFast( unsigned aIndex );


            /// Returns a pointer to the ParamArray (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            const ParamArray* ArrayAt( unsigned aIndex ) const;
            /// Returns a pointer to the ParamArray (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            ParamArray* ArrayAt( unsigned aIndex );
            /// Returns a pointer to the ParamArray (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            const ParamArray* ArrayAtFast( unsigned aIndex ) const;
            /// Returns a pointer to the ParamArray (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            ParamArray* ArrayAtFast( unsigned aIndex );

            /// Returns a pointer to the ParamNode (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            const ParamNode* NodeAt( unsigned aIndex ) const;
            /// Returns a pointer to the ParamNode (dynamic-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present or type is wrong.
            ParamNode* NodeAt( unsigned aIndex );
            /// Returns a pointer to the ParamNode (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            const ParamNode* NodeAtFast( unsigned aIndex ) const;
            /// Returns a pointer to the ParamNode (static-ly cast) corresponding to aName.
            /// Returns NULL if aName is not present.
            ParamNode* NodeAtFast( unsigned aIndex );

            /// Returns a reference to the Param At aIndex.
            /// Behavior is undefined if aIndex is out-of-range.
            const Param& operator[]( unsigned aIndex ) const;
            /// Returns a reference to the Param At aIndex.
            /// Behavior is undefined if aIndex is out-of-range.
            Param& operator[]( unsigned aIndex );

            const Param* Front() const;
            Param* Front();

            const Param* Back() const;
            Param* Back();

            // assign a copy of aValue to the array at aIndex
            void Assign( unsigned aIndex, const Param& aValue );
            // directly assign aValue_ptr to the array at aIndex
            void Assign( unsigned aIndex, Param* aValue_ptr );

            void PushBack( const Param& aValue );
            void PushBack( Param* aValue_ptr );

            void PushFront( const Param& aValue );
            void PushFront( Param* aValue_ptr );

            void Erase( unsigned aIndex );
            Param* Remove( unsigned aIndex );

            iterator Begin();
            const_iterator Begin() const;

            iterator End();
            const_iterator End() const;

            reverse_iterator RBegin();
            const_reverse_iterator RBegin() const;

            reverse_iterator REnd();
            const_reverse_iterator REnd() const;

            virtual std::string ToString() const;

        protected:
            contents fContents;
    };

    template< typename XValType >
    XValType ParamArray::GetValue( unsigned aIndex ) const
    {
        const ParamValue* value = ValueAt( aIndex );
        if( value == NULL ) throw Exception() << "No Value At <" << aIndex << "> is present At this node";
        return value->Get< XValType >();
    }

    template< typename XValType >
    XValType ParamArray::GetValue( unsigned aIndex, XValType aDefault ) const
    {
        const ParamValue* value = ValueAt( aIndex );
        if( value == NULL ) return aDefault;
        return value->Get< XValType >();
    }



    class ParamNode : public Param
    {
        public:
            typedef std::map< std::string, Param* > contents;
            typedef contents::iterator iterator;
            typedef contents::const_iterator const_iterator;
            typedef contents::value_type contents_type;

            ParamNode();
            ParamNode( const ParamNode& orig );
            virtual ~ParamNode();

            virtual Param* Clone() const;

            virtual bool IsNull() const;
            virtual bool IsNode() const;

            bool Has( const std::string& aName ) const;
            unsigned Count( const std::string& aName ) const;

            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Throws a Exception if aName is not present or is not of type ParamValue
            const std::string& GetValue( const std::string& aName ) const;
            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Throws a Exception if aName is not present or is not of type ParamValue
            template< typename XValType >
            XValType GetValue( const std::string& aName ) const;

            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Returns a_default if aName is not present or is not of type ParamValue
            const std::string& GetValue( const std::string& aName, const std::string& aDefault ) const;
            const std::string& GetValue( const std::string& aName, const char* aDefault ) const;
            /// Returns the result of ParamValue::Get if aName is present and is of type ParamValue
            /// Returns a_default if aName is not present or is not of type ParamValue
            template< typename XValType >
            XValType GetValue( const std::string& aName, XValType aDefault ) const;

            /// Returns a pointer to the Param corresponding to aName.
            /// Returns NULL if aName is not present.
            const Param* At( const std::string& aName ) const;
            /// Returns a pointer to the Param corresponding to aName.
            /// Returns NULL if aName is not present.
            Param* At( const std::string& aName );

            const ParamValue* ValueAt( const std::string& aName ) const;
            ParamValue* ValueAt( const std::string& aName );
            const ParamValue* ValueAtFast( const std::string& aName ) const;
            ParamValue* ValueAtFast( const std::string& aName );

            const ParamArray* ArrayAt( const std::string& aName ) const;
            ParamArray* ArrayAt( const std::string& aName );
            const ParamArray* ArrayAtFast( const std::string& aName ) const;
            ParamArray* ArrayAtFast( const std::string& aName );

            const ParamNode* NodeAt( const std::string& aName ) const;
            ParamNode* NodeAt( const std::string& aName );
            const ParamNode* NodeAtFast( const std::string& aName ) const;
            ParamNode* NodeAtFast( const std::string& aName );

            /// Returns a reference to the Param corresponding to aName.
            /// Throws an Exception if aName is not present.
            const Param& operator[]( const std::string& aName ) const;
            /// Returns a reference to the Param corresponding to aName.
            /// Adds a new Value if aName is not present.
            Param& operator[]( const std::string& aName );

            /// creates a copy of aValue
            bool Add( const std::string& aName, const Param& aValue );
            /// directly adds (without copying) aValuePtr
            bool Add( const std::string& aName, Param* aValuePtr );

            /// creates a copy of aValue
            void Replace( const std::string& aName, const Param& aValue );
            /// directly adds (without copying) aValuePtr
            void Replace( const std::string& aName, Param* aValuePtr );

            /// Merges the contents of a_object into this object.
            /// If names in the contents of aObject exist in this object,
            /// the values in this object corresponding to the matching names will be replaced.
            void Merge( const ParamNode& aObject );

            void Erase( const std::string& aName );
            Param* Remove( const std::string& aName );

            iterator Begin();
            const_iterator Begin() const;

            iterator End();
            const_iterator End() const;

            virtual std::string ToString() const;

        protected:
            contents fContents;

    };

    template< typename XValType >
    XValType ParamNode::GetValue( const std::string& aName ) const
    {
        const ParamValue* value = ValueAt( aName );
        if( value == NULL ) throw Exception() << "No Value with name <" << aName << "> is present At this node";
        return value->Get< XValType >();
    }

    template< typename XValType >
    XValType ParamNode::GetValue( const std::string& aName, XValType aDefault ) const
    {
        const ParamValue* value = ValueAt( aName );
        if( value == NULL ) return aDefault;
        return value->Get< XValType >();
    }



    std::ostream& operator<<(std::ostream& out, const Param& aValue);
    std::ostream& operator<<(std::ostream& out, const ParamValue& value);
    std::ostream& operator<<(std::ostream& out, const ParamArray& value);
    std::ostream& operator<<(std::ostream& out, const ParamNode& value);



    //***************************************
    //************** INPUT ******************
    //***************************************

    class ParamInputJSON
    {
        public:
            ParamInputJSON();
            virtual ~ParamInputJSON();

            static ParamNode* ReadFile( const std::string& aFilename );
            static ParamNode* ReadString( const std::string& aJSONStr );
            static ParamNode* ReadDocument( const rapidjson::Document& aDocument );
            static Param* ReadValue( const rapidjson::Value& aValue );
    };

} /* namespace locust */

#endif /* LMCPARAM_HH_ */
