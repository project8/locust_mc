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

            virtual Param* clone() const;

            virtual bool is_null() const;
            virtual bool is_value() const;
            virtual bool is_array() const;
            virtual bool is_node() const;

            ParamValue& as_value();
            ParamArray& as_array();
            ParamNode& as_node();

            const ParamValue& as_value() const;
            const ParamArray& as_array() const;
            const ParamNode& as_node() const;

            virtual std::string to_string() const;

            static unsigned s_indent_level;
    };

    class ParamValue : public Param
    {
        public:
            ParamValue();
            template< typename XStreamableType >
            ParamValue( XStreamableType a_streamable );
            ParamValue(const ParamValue& orig);
            virtual ~ParamValue();

            virtual Param* clone() const;

            virtual bool is_null() const;
            virtual bool is_value() const;

            const std::string& get() const;
            template< typename XValType >
            XValType get() const;

            template< typename XStreamableType >
            ParamValue& operator<<( const XStreamableType& a_streamable );

            virtual std::string to_string() const;

        protected:
            std::string f_value;

    };

    template< typename XStreamableType >
    ParamValue::ParamValue( XStreamableType a_streamable ) :
            Param(),
            f_value()
    {
        (*this) << a_streamable;
    }

    template< typename XValType >
    XValType ParamValue::get() const
    {
        XValType t_return;
        stringstream t_buffer;
        t_buffer << f_value;
        t_buffer >> t_return;
        return t_return;
    }

    template< typename XStreamableType >
    ParamValue& ParamValue::operator<<( const XStreamableType& a_streamable )
    {
        std::stringstream t_buffer;
        t_buffer << a_streamable;
        f_value = t_buffer.str();
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

            virtual Param* clone() const;

            virtual bool is_null() const;
            virtual bool is_array() const;

            unsigned size() const;
            bool empty() const;

            /// sets the size of the array
            /// if smaller than the current size, extra elements are deleted
            void resize( unsigned a_size );

            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Throws an exception if a_name is not present or is not of type ParamValue
            const std::string& get_value( unsigned a_index ) const;
            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Throws an exception if a_name is not present or is not of type ParamValue
            template< typename XValType >
            XValType get_value( unsigned a_index ) const;

            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Returns a_default if a_name is not present or is not of type ParamValue
            const std::string& get_value( unsigned a_index, const std::string& a_default ) const;
            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Returns a_default if a_name is not present or is not of type ParamValue
            template< typename XValType >
            XValType get_value( unsigned a_index, XValType a_default ) const;

            /// Returns a pointer to the Param corresponding to a_name.
            /// Returns NULL if a_name is not present.
            const Param* at( unsigned a_index ) const;
            /// Returns a pointer to the Param corresponding to a_name.
            /// Returns NULL if a_name is not present.
            Param* at( unsigned a_index );

            /// Returns a pointer to the ParamValue (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            const ParamValue* value_at( unsigned a_index ) const;
            /// Returns a pointer to the ParamValue (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            ParamValue* value_at( unsigned a_index );

            /// Returns a pointer to the ParamArray (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            const ParamArray* array_at( unsigned a_index ) const;
            /// Returns a pointer to the ParamArray (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            ParamArray* array_at( unsigned a_index );

            /// Returns a pointer to the ParamNode (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            const ParamNode* node_at( unsigned a_index ) const;
            /// Returns a pointer to the ParamNode (static-ly cast) corresponding to a_name.
            /// Returns NULL if a_name is not present.
            ParamNode* node_at( unsigned a_index );

            /// Returns a reference to the Param at a_index.
            /// Behavior is undefined if a_index is out-of-range.
            const Param& operator[]( unsigned a_index ) const;
            /// Returns a reference to the Param at a_index.
            /// Behavior is undefined if a_index is out-of-range.
            Param& operator[]( unsigned a_index );

            const Param* front() const;
            Param* front();

            const Param* back() const;
            Param* back();

            // assign a copy of a_value to the array at a_index
            void assign( unsigned a_index, const Param& a_value );
            // directly assign a_value_ptr to the array at a_index
            void assign( unsigned a_index, Param* a_value_ptr );

            void push_back( const Param& a_value );
            void push_back( Param* a_value_ptr );

            void push_front( const Param& a_value );
            void push_front( Param* a_value_ptr );

            void erase( unsigned a_index );
            Param* remove( unsigned a_index );

            iterator begin();
            const_iterator begin() const;

            iterator end();
            const_iterator end() const;

            reverse_iterator rbegin();
            const_reverse_iterator rbegin() const;

            reverse_iterator rend();
            const_reverse_iterator rend() const;

            virtual std::string to_string() const;

        protected:
            contents f_contents;
    };

    template< typename XValType >
    XValType ParamArray::get_value( unsigned a_index ) const
    {
        const ParamValue* value = value_at( a_index );
        if( value == NULL ) throw Exception() << "No value at <" << a_index << "> is present at this node";
        return value->get< XValType >();
    }

    template< typename XValType >
    XValType ParamArray::get_value( unsigned a_index, XValType a_default ) const
    {
        const ParamValue* value = value_at( a_index );
        if( value == NULL ) return a_default;
        return value->get< XValType >();
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

            virtual Param* clone() const;

            virtual bool is_null() const;
            virtual bool is_node() const;

            bool has( const std::string& a_name ) const;
            unsigned count( const std::string& a_name ) const;

            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Throws an exception if a_name is not present or is not of type ParamValue
            const std::string& get_value( const std::string& a_name ) const;
            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Throws an exception if a_name is not present or is not of type ParamValue
            template< typename XValType >
            XValType get_value( const std::string& a_name ) const;

            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Returns a_default if a_name is not present or is not of type ParamValue
            const std::string& get_value( const std::string& a_name, const std::string& a_default ) const;
            /// Returns the result of ParamValue::get if a_name is present and is of type ParamValue
            /// Returns a_default if a_name is not present or is not of type ParamValue
            template< typename XValType >
            XValType get_value( const std::string& a_name, XValType a_default ) const;

            /// Returns a pointer to the Param corresponding to a_name.
            /// Returns NULL if a_name is not present.
            const Param* at( const std::string& a_name ) const;
            /// Returns a pointer to the Param corresponding to a_name.
            /// Returns NULL if a_name is not present.
            Param* at( const std::string& a_name );

            const ParamValue* value_at( const std::string& a_name ) const;
            ParamValue* value_at( const std::string& a_name );

            const ParamArray* array_at( const std::string& a_name ) const;
            ParamArray* array_at( const std::string& a_name );

            const ParamNode* node_at( const std::string& a_name ) const;
            ParamNode* node_at( const std::string& a_name );

            /// Returns a reference to the Param corresponding to a_name.
            /// Throws an exception if a_name is not present.
            const Param& operator[]( const std::string& a_name ) const;
            /// Returns a reference to the Param corresponding to a_name.
            /// Adds a new value if a_name is not present.
            Param& operator[]( const std::string& a_name );

            /// creates a copy of a_value
            bool add( const std::string& a_name, const Param& a_value );
            /// directly adds (without copying) a_value_ptr
            bool add( const std::string& a_name, Param* a_value_ptr );

            /// creates a copy of a_value
            void replace( const std::string& a_name, const Param& a_value );
            /// directly adds (without copying) a_value_ptr
            void replace( const std::string& a_name, Param* a_value_ptr );

            /// Merges the contents of a_object into this object.
            /// If names in the contents of a_object exist in this object,
            /// the values in this object corresponding to the matching names will be replaced.
            void merge( const ParamNode* a_object );

            void erase( const std::string& a_name );
            Param* remove( const std::string& a_name );

            iterator begin();
            const_iterator begin() const;

            iterator end();
            const_iterator end() const;

            virtual std::string to_string() const;

        protected:
            contents f_contents;

    };

    template< typename XValType >
    XValType ParamNode::get_value( const std::string& a_name ) const
    {
        const ParamValue* value = value_at( a_name );
        if( value == NULL ) throw Exception() << "No value with name <" << a_name << "> is present at this node";
        return value->get< XValType >();
    }

    template< typename XValType >
    XValType ParamNode::get_value( const std::string& a_name, XValType a_default ) const
    {
        const ParamValue* value = value_at( a_name );
        if( value == NULL ) return a_default;
        return value->get< XValType >();
    }




    std::ostream& operator<<(std::ostream& out, const Param& value);
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

            static ParamNode* read_file( const std::string& a_filename );
            static ParamNode* read_string( const std::string& a_json_str );
            static ParamNode* read_document( const rapidjson::Document& a_document );
            static Param* read_value( const rapidjson::Value& a_value );
    };

} /* namespace locust */

#endif /* LMCPARAM_HH_ */
