/*
 * LMCParam.cc
 *
 *  Created on: Jan 14, 2014
 *      Author: nsoblath
 */

#include "LMCParam.hh"

#include "LMCLogger.hh"

#include <sstream>
using std::string;
using std::stringstream;


// for the config maker
#include "filestream.h"

#include <cstdio>
using std::string;



namespace locust
{
    LMCLOGGER( lmclog, "Param" );

    unsigned Param::s_indent_level = 0;

    Param::Param()
    {
    }

    Param::Param( const Param& )
    {
    }

    Param::~Param()
    {
    }

    Param* Param::clone() const
    {
        return new Param( *this );
    }

    bool Param::is_null() const
    {
        return true;
    }

    bool Param::is_value() const
    {
        return false;
    }

    bool Param::is_array() const
    {
        return false;
    }

    bool Param::is_node() const
    {
        return false;
    }

    ParamValue& Param::as_value()
    {
        ParamValue* t_cast_ptr = static_cast< ParamValue* >( this );
        return *t_cast_ptr;
    }

    ParamArray& Param::as_array()
    {
        ParamArray* t_cast_ptr = static_cast< ParamArray* >( this );
        return *t_cast_ptr;
    }

    ParamNode& Param::as_node()
    {
        ParamNode* t_cast_ptr = static_cast< ParamNode* >( this );
        return *t_cast_ptr;
    }

    const ParamValue& Param::as_value() const
    {
        const ParamValue* t_cast_ptr = dynamic_cast< const ParamValue* >( this );
        return *t_cast_ptr;
    }

    const ParamArray& Param::as_array() const
    {
        const ParamArray* t_cast_ptr = dynamic_cast< const ParamArray* >( this );
        return *t_cast_ptr;
    }

    const ParamNode& Param::as_node() const
    {
        const ParamNode* t_cast_ptr = dynamic_cast< const ParamNode* >( this );
        return *t_cast_ptr;
    }

    std::string Param::to_string() const
    {
        return string();
    }

    //************************************
    //***********  DATA  *****************
    //************************************

    ParamValue::ParamValue() :
            Param(),
            f_value_str(),
            f_value_str_buffer(),
            f_value_buffer()
    {
    }

    ParamValue::ParamValue( const ParamValue& orig ) :
            Param( orig ),
            f_value_str()
    {
        f_value_str << orig.f_value_str.str();
    }

    ParamValue::~ParamValue()
    {
    }

    Param* ParamValue::clone() const
    {
        return new ParamValue( *this );
    }

    bool ParamValue::is_null() const
    {
        return false;
    }

    bool ParamValue::is_value() const
    {
        return true;
    }

    const string& ParamValue::get() const
    {
        f_value_buffer = f_value_str.str();
        return f_value_buffer;
    }

    std::string ParamValue::to_string() const
    {
        return string(f_value_str.str());
    }

    //************************************
    //***********  ARRAY  ****************
    //************************************

    ParamArray::ParamArray() :
            Param(),
            f_contents()
    {
    }

    ParamArray::ParamArray( const ParamArray& orig ) :
            Param( orig ),
            f_contents( orig.f_contents.size() )
    {
        for( unsigned ind = 0; ind < f_contents.size(); ++ind )
        {
            this->assign( ind, orig[ ind ].clone() );
        }
    }

    ParamArray::~ParamArray()
    {
        for( unsigned ind = 0; ind < f_contents.size(); ++ind )
        {
            delete f_contents[ ind ];
        }
    }

    Param* ParamArray::clone() const
    {
        return new ParamArray( *this );
    }

    bool ParamArray::is_null() const
    {
        return false;
    }

    bool ParamArray::is_array() const
    {
        return true;
    }

    unsigned ParamArray::size() const
    {
        return f_contents.size();
    }
    bool ParamArray::empty() const
    {
        return f_contents.empty();
    }

    void ParamArray::resize( unsigned a_size )
    {
        unsigned curr_size = f_contents.size();
        for( unsigned ind = a_size; ind < curr_size; ++ind )
        {
            delete f_contents[ ind ];
        }
        f_contents.resize( a_size );
        return;
    }

    const std::string& ParamArray::get_value( unsigned a_index ) const
    {
        const ParamValue* value = value_at( a_index );
        if( value == NULL ) throw Exception() << "No value at <" << a_index << "> is present at this node";
        return value->get();
    }

    const std::string& ParamArray::get_value( unsigned a_index, const std::string& a_default ) const
    {
        const ParamValue* value = value_at( a_index );
        if( value == NULL ) return a_default;
        return value->get();
    }

    const Param* ParamArray::at( unsigned a_index ) const
    {
        if( a_index >= f_contents.size() ) return NULL;
        return f_contents[ a_index ];
    }
    Param* ParamArray::at( unsigned a_index )
    {
        if( a_index >= f_contents.size() ) return NULL;
        return f_contents[ a_index ];
    }

    const ParamValue* ParamArray::value_at( unsigned a_index ) const
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_value();
    }
    ParamValue* ParamArray::value_at( unsigned a_index )
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_value();
    }

    const ParamArray* ParamArray::array_at( unsigned a_index ) const
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_array();
    }
    ParamArray* ParamArray::array_at( unsigned a_index )
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_array();
    }

    const ParamNode* ParamArray::node_at( unsigned a_index ) const
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_node();
    }
    ParamNode* ParamArray::node_at( unsigned a_index )
    {
        if( a_index >= f_contents.size() ) return NULL;
        return &f_contents[ a_index ]->as_node();
    }

    const Param& ParamArray::operator[]( unsigned a_index ) const
    {
        return *f_contents[ a_index ];
    }
    Param& ParamArray::operator[]( unsigned a_index )
    {
        return *f_contents[ a_index ];
    }

    const Param* ParamArray::front() const
    {
        return f_contents.front();
    }
    Param* ParamArray::front()
    {
        return f_contents.front();
    }

    const Param* ParamArray::back() const
    {
        return f_contents.back();
    }
    Param* ParamArray::back()
    {
        return f_contents.back();
    }

    // assign a copy of a_value to the array at a_index
    void ParamArray::assign( unsigned a_index, const Param& a_value )
    {
        erase( a_index );
        f_contents[ a_index ] = a_value.clone();
        return;
    }
    // directly assign a_value_ptr to the array at a_index
    void ParamArray::assign( unsigned a_index, Param* a_value_ptr )
    {
        erase( a_index );
        f_contents[ a_index ] = a_value_ptr;
        return;
    }

    void ParamArray::push_back( const Param& a_value )
    {
        f_contents.push_back( a_value.clone() );
        return;
    }
    void ParamArray::push_back( Param* a_value_ptr )
    {
        f_contents.push_back( a_value_ptr );
        return;
    }

    void ParamArray::push_front( const Param& a_value )
    {
        f_contents.push_front( a_value.clone() );
        return;
    }
    void ParamArray::push_front( Param* a_value_ptr )
    {
        f_contents.push_front( a_value_ptr );
        return;
    }

    void ParamArray::erase( unsigned a_index )
    {
        delete f_contents[ a_index ];
        return;
    }
    Param* ParamArray::remove( unsigned a_index )
    {
        Param* t_current = f_contents[ a_index ];
        f_contents[ a_index ] = NULL;
        return t_current;
    }

    ParamArray::iterator ParamArray::begin()
    {
        return f_contents.begin();
    }
    ParamArray::const_iterator ParamArray::begin() const
    {
        return f_contents.begin();
    }

    ParamArray::iterator ParamArray::end()
    {
        return f_contents.end();
    }
    ParamArray::const_iterator ParamArray::end() const
    {
        return f_contents.end();
    }

    ParamArray::reverse_iterator ParamArray::rbegin()
    {
        return f_contents.rbegin();
    }
    ParamArray::const_reverse_iterator ParamArray::rbegin() const
    {
        return f_contents.rbegin();
    }

    ParamArray::reverse_iterator ParamArray::rend()
    {
        return f_contents.rend();
    }
    ParamArray::const_reverse_iterator ParamArray::rend() const
    {
        return f_contents.rend();
    }

    std::string ParamArray::to_string() const
    {
        stringstream out;
        string indentation;
        for ( unsigned i=0; i<Param::s_indent_level; ++i )
            indentation += "    ";
        out << '\n' << indentation << "[\n";
        Param::s_indent_level++;
        for( const_iterator it = begin(); it != end(); ++it )
        {
            out << indentation << "    " << **it << '\n';
        }
        Param::s_indent_level--;
        out << indentation << "]\n";
        return out.str();
    }


    //************************************
    //***********  NODE  *****************
    //************************************

    ParamNode::ParamNode() :
            Param(),
            f_contents()
    {
    }

    ParamNode::ParamNode( const ParamNode& orig ) :
            Param( orig ),
            f_contents()
    {
        for( const_iterator it = orig.f_contents.begin(); it != orig.f_contents.end(); ++it )
        {
            this->replace( it->first, *it->second );
        }
    }

    ParamNode::~ParamNode()
    {
        for( iterator it = f_contents.begin(); it != f_contents.end(); ++it )
        {
            delete it->second;
        }
    }

    Param* ParamNode::clone() const
    {
        return new ParamNode( *this );
    }

    bool ParamNode::is_null() const
    {
        return false;
    }

    bool ParamNode::is_node() const
    {
        return true;
    }

    bool ParamNode::has( const std::string& a_name ) const
    {
        return f_contents.count( a_name ) > 0;
    }

    unsigned ParamNode::count( const std::string& a_name ) const
    {
        return f_contents.count( a_name );
    }

    const std::string& ParamNode::get_value( const std::string& a_name ) const
    {
        const ParamValue* value = value_at( a_name );
        if( value == NULL ) throw Exception() << "No value with name <" << a_name << "> is present at this node";
        return value->get();
    }

    const std::string& ParamNode::get_value( const std::string& a_name, const std::string& a_default ) const
    {
        const ParamValue* value = value_at( a_name );
        if( value == NULL ) return a_default;
        return value->get();
    }

    const Param* ParamNode::at( const std::string& a_name ) const
    {
        const_iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return it->second;
    }

    Param* ParamNode::at( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return it->second;
    }

    const ParamValue* ParamNode::value_at( const std::string& a_name ) const
    {
        const_iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_value();
    }

    ParamValue* ParamNode::value_at( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_value();
    }

    const ParamArray* ParamNode::array_at( const std::string& a_name ) const
    {
        const_iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_array();
    }

    ParamArray* ParamNode::array_at( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_array();
    }

    const ParamNode* ParamNode::node_at( const std::string& a_name ) const
    {

        const_iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_node();
    }

    ParamNode* ParamNode::node_at( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            return NULL;
        }
        return &it->second->as_node();
    }

    const Param& ParamNode::operator[]( const std::string& a_name ) const
    {
        const_iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            throw Exception() << "No value present corresponding to name <" << a_name << ">\n";
        }
        return *(it->second);
    }

    Param& ParamNode::operator[]( const std::string& a_name )
    {
        return *f_contents[ a_name ];
    }

    bool ParamNode::add( const std::string& a_name, const Param& a_value )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            f_contents.insert( contents_type( a_name, a_value.clone() ) );
            return true;
        }
        return false;
    }

    bool ParamNode::add( const std::string& a_name, Param* a_value )
    {
        iterator it = f_contents.find( a_name );
        if( it == f_contents.end() )
        {
            f_contents.insert( contents_type( a_name, a_value ) );
            return true;
        }
        return false;
    }

    void ParamNode::replace( const std::string& a_name, const Param& a_value )
    {
        erase( a_name );
        f_contents[ a_name ] = a_value.clone();
        return;
    }

    void ParamNode::replace( const std::string& a_name, Param* a_value )
    {
        erase( a_name );
        f_contents[ a_name ] = a_value;
        return;
    }

    void ParamNode::merge( const ParamNode* a_object )
    {
        for( const_iterator it = a_object->f_contents.begin(); it != a_object->f_contents.end(); ++it )
        {
            this->replace( it->first, *it->second );
        }
    }

    void ParamNode::erase( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it != f_contents.end() )
        {
            delete it->second;
            f_contents.erase( it );
        }
        return;
    }

    Param* ParamNode::remove( const std::string& a_name )
    {
        iterator it = f_contents.find( a_name );
        if( it != f_contents.end() )
        {
            Param* removed = it->second;
            f_contents.erase( it );
            return removed;
        }
        return NULL;
    }

    ParamNode::iterator ParamNode::begin()
    {
        return f_contents.begin();
    }

    ParamNode::const_iterator ParamNode::begin() const
    {
        return f_contents.begin();
    }

    ParamNode::iterator ParamNode::end()
    {
        return f_contents.end();
    }

    ParamNode::const_iterator ParamNode::end() const
    {
        return f_contents.end();
    }

    std::string ParamNode::to_string() const
    {
        stringstream out;
        string indentation;
        for ( unsigned i=0; i<Param::s_indent_level; ++i )
            indentation += "    ";
        out << '\n' << indentation << "{\n";
        Param::s_indent_level++;
        for( const_iterator it = begin(); it != end(); ++it )
        {
            out << indentation << "    " << it->first << " : " << *(it->second) << '\n';
        }
        Param::s_indent_level--;
        out << indentation << "}\n";
        return out.str();
    }




    std::ostream& operator<<(std::ostream& out, const Param& a_value)
    {
        return out << a_value.to_string();
    }


    std::ostream& operator<<(std::ostream& out, const ParamValue& a_value)
    {
        return out << a_value.to_string();
    }


    std::ostream& operator<<(std::ostream& out, const ParamArray& a_value)
    {
        return out << a_value.to_string();
    }


    std::ostream& operator<<(std::ostream& out, const ParamNode& a_value)
    {
        return out << a_value.to_string();
    }





    ParamInputJSON::ParamInputJSON()
    {
    }
    ParamInputJSON::~ParamInputJSON()
    {
    }

    ParamNode* ParamInputJSON::read_file( const std::string& a_filename )
    {
        FILE* t_config_file = fopen( a_filename.c_str(), "r" );
        if( t_config_file == NULL )
        {
            LMCERROR( lmclog, "file <" << a_filename << "> did not open" );
            return NULL;
        }
        rapidjson::FileStream t_file_stream( t_config_file );

        rapidjson::Document t_config_doc;
        if( t_config_doc.ParseStream<0>( t_file_stream ).HasParseError() )
        {
            unsigned errorPos = t_config_doc.GetErrorOffset();
            rewind( t_config_file );
            unsigned iChar, newlineCount = 1, lastNewlinePos = 0;
            int thisChar;
            for( iChar = 0; iChar != errorPos; ++iChar )
            {
                thisChar = fgetc( t_config_file );
                if( thisChar == EOF )
                {
                    break;
                }
                if( thisChar == '\n' || thisChar == '\r' )
                {
                    newlineCount++;
                    lastNewlinePos = iChar + 1;
                }
            }
            if( iChar == errorPos )
            {
                LMCERROR( lmclog, "error parsing config file :\n" <<
                        '\t' << t_config_doc.GetParseError() << '\n' <<
                        "\tThe error was reported at line " << newlineCount << ", character " << errorPos - lastNewlinePos );
            }
            else
            {
                LMCERROR( lmclog, "error parsing config file :\n" <<
                        '\t' << t_config_doc.GetParseError() <<
                        "\tend of file reached before error location was found" );
            }
            fclose( t_config_file );
            return NULL;
        }
        fclose( t_config_file );

        return ParamInputJSON::read_document( t_config_doc );
    }

    ParamNode* ParamInputJSON::read_string( const std::string& a_json_string )
    {
        rapidjson::Document t_config_doc;
        if( t_config_doc.Parse<0>( a_json_string.c_str() ).HasParseError() )
        {
            LMCERROR( lmclog, "error parsing string:\n" << t_config_doc.GetParseError() );
            return NULL;
        }
        return ParamInputJSON::read_document( t_config_doc );
    }

    ParamNode* ParamInputJSON::read_document( const rapidjson::Document& a_doc )
    {
        ParamNode* t_config = new ParamNode();
        for( rapidjson::Value::ConstMemberIterator jsonIt = a_doc.MemberBegin();
                jsonIt != a_doc.MemberEnd();
                ++jsonIt)
        {
            t_config->replace( jsonIt->name.GetString(), ParamInputJSON::read_value( jsonIt->value ) );
        }
        return t_config;
    }

    Param* ParamInputJSON::read_value( const rapidjson::Value& a_value )
    {
        if( a_value.IsNull() )
        {
            return new Param();
        }
        if( a_value.IsObject() )
        {
            ParamNode* t_config_object = new ParamNode();
            for( rapidjson::Value::ConstMemberIterator jsonIt = a_value.MemberBegin();
                    jsonIt != a_value.MemberEnd();
                    ++jsonIt)
            {
                t_config_object->replace( jsonIt->name.GetString(), ParamInputJSON::read_value( jsonIt->value ) );
            }
            return t_config_object;
        }
        if( a_value.IsArray() )
        {
            ParamArray* t_config_array = new ParamArray();
            for( rapidjson::Value::ConstValueIterator jsonIt = a_value.Begin();
                    jsonIt != a_value.End();
                    ++jsonIt)
            {
                t_config_array->push_back( ParamInputJSON::read_value( *jsonIt ) );
            }
            return t_config_array;
        }
        if( a_value.IsString() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetString();
            return t_config_value;
        }
        if( a_value.IsBool() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetBool();
            return t_config_value;
        }
        if( a_value.IsInt() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetInt();
            return t_config_value;
        }
        if( a_value.IsUint() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetUint();
            return t_config_value;
        }
        if( a_value.IsInt64() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetInt64();
            return t_config_value;
        }
        if( a_value.IsUint64() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetUint64();
            return t_config_value;
        }
        if( a_value.IsDouble() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << a_value.GetDouble();
            return t_config_value;
        }
        LMCWARN( lmclog, "unknown type; returning null value" );
        return new Param();
    }

} /* namespace locust */
