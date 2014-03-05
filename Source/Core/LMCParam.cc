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
    LMCLOGGER( plog, "Param" );

    unsigned Param::sIndentLevel = 0;

    Param::Param()
    {
    }

    Param::Param( const Param& )
    {
    }

    Param::~Param()
    {
    }

    Param* Param::Clone() const
    {
        return new Param( *this );
    }

    bool Param::IsNull() const
    {
        return true;
    }

    bool Param::IsValue() const
    {
        return false;
    }

    bool Param::IsArray() const
    {
        return false;
    }

    bool Param::IsNode() const
    {
        return false;
    }

    ParamValue& Param::AsValue()
    {
        ParamValue* t_cast_ptr = static_cast< ParamValue* >( this );
        return *t_cast_ptr;
    }

    ParamArray& Param::AsArray()
    {
        ParamArray* t_cast_ptr = static_cast< ParamArray* >( this );
        return *t_cast_ptr;
    }

    ParamNode& Param::AsNode()
    {
        ParamNode* t_cast_ptr = static_cast< ParamNode* >( this );
        return *t_cast_ptr;
    }

    const ParamValue& Param::AsValue() const
    {
        const ParamValue* t_cast_ptr = static_cast< const ParamValue* >( this );
        return *t_cast_ptr;
    }

    const ParamArray& Param::AsArray() const
    {
        const ParamArray* t_cast_ptr = static_cast< const ParamArray* >( this );
        return *t_cast_ptr;
    }

    const ParamNode& Param::AsNode() const
    {
        const ParamNode* t_cast_ptr = static_cast< const ParamNode* >( this );
        return *t_cast_ptr;
    }

    const ParamValue& Param::operator()() const
    {
        return AsValue();
    }

    ParamValue& Param::operator()()
    {
        return AsValue();
    }

    const Param& Param::operator[]( unsigned aIndex ) const
    {
        return AsArray()[aIndex];
    }

    Param& Param::operator[]( unsigned aIndex )
    {
        return AsArray()[aIndex];
    }

    const Param& Param::operator[]( const std::string& aName ) const
    {
        return AsNode()[aName];
    }

    Param& Param::operator[]( const std::string& aName )
    {
        return AsNode()[aName];
    }

    std::string Param::ToString() const
    {
        return string();
    }

    //************************************
    //***********  DATA  *****************
    //************************************

    ParamValue::ParamValue() :
            Param(),
            fValue()
    {
    }

    ParamValue::ParamValue( const ParamValue& orig ) :
            Param( orig ),
            fValue( orig.fValue )
    {
    }

    ParamValue::~ParamValue()
    {
    }

    Param* ParamValue::Clone() const
    {
        return new ParamValue( *this );
    }

    bool ParamValue::IsNull() const
    {
        return false;
    }

    bool ParamValue::IsValue() const
    {
        return true;
    }

    const string& ParamValue::Get() const
    {
         return fValue;
    }

    std::string ParamValue::ToString() const
    {
        return string( fValue );
    }

    //************************************
    //***********  ARRAY  ****************
    //************************************

    ParamArray::ParamArray() :
            Param(),
            fContents()
    {
    }

    ParamArray::ParamArray( const ParamArray& orig ) :
            Param( orig ),
            fContents( orig.fContents.size() )
    {
        for( unsigned ind = 0; ind < fContents.size(); ++ind )
        {
            this->Assign( ind, orig[ ind ].Clone() );
        }
    }

    ParamArray::~ParamArray()
    {
        for( unsigned ind = 0; ind < fContents.size(); ++ind )
        {
            delete fContents[ ind ];
        }
    }

    Param* ParamArray::Clone() const
    {
        return new ParamArray( *this );
    }

    bool ParamArray::IsNull() const
    {
        return false;
    }

    bool ParamArray::IsArray() const
    {
        return true;
    }

    unsigned ParamArray::Size() const
    {
        return fContents.size();
    }
    bool ParamArray::Empty() const
    {
        return fContents.empty();
    }

    void ParamArray::Resize( unsigned a_size )
    {
        unsigned curr_size = fContents.size();
        for( unsigned ind = a_size; ind < curr_size; ++ind )
        {
            delete fContents[ ind ];
        }
        fContents.resize( a_size );
        return;
    }

    const std::string& ParamArray::GetValue( unsigned aIndex ) const
    {
        const ParamValue* value = ValueAt( aIndex );
        if( value == NULL ) throw Exception() << "No value at <" << aIndex << "> is present at this node";
        return value->Get();
    }

    const std::string& ParamArray::GetValue( unsigned aIndex, const std::string& aDefault ) const
    {
        const ParamValue* value = ValueAt( aIndex );
        if( value == NULL ) return aDefault;
        return value->Get();
    }

    const std::string& ParamArray::GetValue( unsigned aIndex, const char* aDefault ) const
    {
        return GetValue( aIndex, string( aDefault ) );
    }

    const Param* ParamArray::At( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return fContents[ aIndex ];
    }
    Param* ParamArray::At( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return fContents[ aIndex ];
    }

    const ParamValue* ParamArray::ValueAt( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamValue* >( fContents[ aIndex ] );
    }
    ParamValue* ParamArray::ValueAt( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamValue* >( fContents[ aIndex ] );
    }

    const ParamValue* ParamArray::ValueAtFast( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamValue* >( fContents[ aIndex ] );
    }
    ParamValue* ParamArray::ValueAtFast( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamValue* >( fContents[ aIndex ] );
    }

    const ParamArray* ParamArray::ArrayAt( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamArray* >( fContents[ aIndex ] );
    }
    ParamArray* ParamArray::ArrayAt( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamArray* >( fContents[ aIndex ] );
    }

    const ParamArray* ParamArray::ArrayAtFast( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamArray* >( fContents[ aIndex ] );
    }
    ParamArray* ParamArray::ArrayAtFast( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamArray* >( fContents[ aIndex ] );
    }

    const ParamNode* ParamArray::NodeAt( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamNode* >( fContents[ aIndex ] );
    }
    ParamNode* ParamArray::NodeAt( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return dynamic_cast< ParamNode* >( fContents[ aIndex ] );
    }

    const ParamNode* ParamArray::NodeAtFast( unsigned aIndex ) const
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamNode* >( fContents[ aIndex ] );
    }
    ParamNode* ParamArray::NodeAtFast( unsigned aIndex )
    {
        if( aIndex >= fContents.size() ) return NULL;
        return static_cast< ParamNode* >( fContents[ aIndex ] );
    }

    const Param& ParamArray::operator[]( unsigned aIndex ) const
    {
        return *fContents[ aIndex ];
    }
    Param& ParamArray::operator[]( unsigned aIndex )
    {
        return *fContents[ aIndex ];
    }

    const Param* ParamArray::Front() const
    {
        return fContents.front();
    }
    Param* ParamArray::Front()
    {
        return fContents.front();
    }

    const Param* ParamArray::Back() const
    {
        return fContents.back();
    }
    Param* ParamArray::Back()
    {
        return fContents.back();
    }

    // assign a copy of aValue to the array at aIndex
    void ParamArray::Assign( unsigned aIndex, const Param& aValue )
    {
        Erase( aIndex );
        fContents[ aIndex ] = aValue.Clone();
        return;
    }
    // directly assign aValue_ptr to the array at aIndex
    void ParamArray::Assign( unsigned aIndex, Param* aValue_ptr )
    {
        Erase( aIndex );
        fContents[ aIndex ] = aValue_ptr;
        return;
    }

    void ParamArray::PushBack( const Param& aValue )
    {
        fContents.push_back( aValue.Clone() );
        return;
    }
    void ParamArray::PushBack( Param* aValue_ptr )
    {
        fContents.push_back( aValue_ptr );
        return;
    }

    void ParamArray::PushFront( const Param& aValue )
    {
        fContents.push_front( aValue.Clone() );
        return;
    }
    void ParamArray::PushFront( Param* aValue_ptr )
    {
        fContents.push_front( aValue_ptr );
        return;
    }

    void ParamArray::Erase( unsigned aIndex )
    {
        delete fContents[ aIndex ];
        return;
    }
    Param* ParamArray::Remove( unsigned aIndex )
    {
        Param* tCurrent = fContents[ aIndex ];
        fContents[ aIndex ] = NULL;
        return tCurrent;
    }

    ParamArray::iterator ParamArray::Begin()
    {
        return fContents.begin();
    }
    ParamArray::const_iterator ParamArray::Begin() const
    {
        return fContents.begin();
    }

    ParamArray::iterator ParamArray::End()
    {
        return fContents.end();
    }
    ParamArray::const_iterator ParamArray::End() const
    {
        return fContents.end();
    }

    ParamArray::reverse_iterator ParamArray::RBegin()
    {
        return fContents.rbegin();
    }
    ParamArray::const_reverse_iterator ParamArray::RBegin() const
    {
        return fContents.rbegin();
    }

    ParamArray::reverse_iterator ParamArray::REnd()
    {
        return fContents.rend();
    }
    ParamArray::const_reverse_iterator ParamArray::REnd() const
    {
        return fContents.rend();
    }

    std::string ParamArray::ToString() const
    {
        stringstream out;
        string indentation;
        for ( unsigned i=0; i<Param::sIndentLevel; ++i )
            indentation += "    ";
        out << '\n' << indentation << "[\n";
        Param::sIndentLevel++;
        for( const_iterator it = Begin(); it != End(); ++it )
        {
            out << indentation << "    " << **it << '\n';
        }
        Param::sIndentLevel--;
        out << indentation << "]";
        return out.str();
    }


    //************************************
    //***********  NODE  *****************
    //************************************

    ParamNode::ParamNode() :
            Param(),
            fContents()
    {
    }

    ParamNode::ParamNode( const ParamNode& orig ) :
            Param( orig ),
            fContents()
    {
        for( const_iterator it = orig.fContents.begin(); it != orig.fContents.end(); ++it )
        {
            this->Replace( it->first, *it->second );
        }
    }

    ParamNode::~ParamNode()
    {
        for( iterator it = fContents.begin(); it != fContents.end(); ++it )
        {
            delete it->second;
        }
    }

    Param* ParamNode::Clone() const
    {
        return new ParamNode( *this );
    }

    bool ParamNode::IsNull() const
    {
        return false;
    }

    bool ParamNode::IsNode() const
    {
        return true;
    }

    bool ParamNode::Has( const std::string& aName ) const
    {
        return fContents.count( aName ) > 0;
    }

    unsigned ParamNode::Count( const std::string& aName ) const
    {
        return fContents.count( aName );
    }

    const std::string& ParamNode::GetValue( const std::string& aName ) const
    {
        const ParamValue* value = ValueAt( aName );
        if( value == NULL ) throw Exception() << "No value with name <" << aName << "> is present at this node";
        return value->Get();
    }

    const std::string& ParamNode::GetValue( const std::string& aName, const std::string& aDefault ) const
    {
        const ParamValue* value = ValueAt( aName );
        if( value == NULL ) return aDefault;
        return value->Get();
    }

    const std::string& ParamNode::GetValue( const std::string& aName, const char* aDefault ) const
    {
        return GetValue( aName, string( aDefault ) );
    }

    const Param* ParamNode::At( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() )
        {
            return NULL;
        }
        return it->second;
    }

    Param* ParamNode::At( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() )
        {
            return NULL;
        }
        return it->second;
    }

    const ParamValue* ParamNode::ValueAt( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamValue* >( it->second );
    }

    ParamValue* ParamNode::ValueAt( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamValue* >( it->second );
    }

    const ParamValue* ParamNode::ValueAtFast( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamValue* >( it->second );
    }

    ParamValue* ParamNode::ValueAtFast( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamValue* >( it->second );
    }

    const ParamArray* ParamNode::ArrayAt( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamArray* >( it->second );
    }

    ParamArray* ParamNode::ArrayAt( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamArray* >( it->second );
    }

    const ParamArray* ParamNode::ArrayAtFast( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamArray* >( it->second );
    }

    ParamArray* ParamNode::ArrayAtFast( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamArray* >( it->second );
    }

    const ParamNode* ParamNode::NodeAt( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamNode* >( it->second );
    }

    ParamNode* ParamNode::NodeAt( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return dynamic_cast< ParamNode* >( it->second );
    }

    const ParamNode* ParamNode::NodeAtFast( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamNode* >( it->second );
    }

    ParamNode* ParamNode::NodeAtFast( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() ) return NULL;
        return static_cast< ParamNode* >( it->second );
    }

    const Param& ParamNode::operator[]( const std::string& aName ) const
    {
        const_iterator it = fContents.find( aName );
        if( it == fContents.end() )
        {
            throw Exception() << "No value present corresponding to name <" << aName << ">\n";
        }
        return *(it->second);
    }

    Param& ParamNode::operator[]( const std::string& aName )
    {
        return *fContents[ aName ];
    }

    bool ParamNode::Add( const std::string& aName, const Param& aValue )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() )
        {
            fContents.insert( contents_type( aName, aValue.Clone() ) );
            return true;
        }
        return false;
    }

    bool ParamNode::Add( const std::string& aName, Param* aValue )
    {
        iterator it = fContents.find( aName );
        if( it == fContents.end() )
        {
            fContents.insert( contents_type( aName, aValue ) );
            return true;
        }
        return false;
    }

    void ParamNode::Replace( const std::string& aName, const Param& aValue )
    {
        Erase( aName );
        fContents[ aName ] = aValue.Clone();
        return;
    }

    void ParamNode::Replace( const std::string& aName, Param* aValue )
    {
        Erase( aName );
        fContents[ aName ] = aValue;
        return;
    }

    void ParamNode::Merge( const ParamNode& aObject )
    {
        for( const_iterator it = aObject.fContents.begin(); it != aObject.fContents.end(); ++it )
        {
            if( this->Has( it->first ) )
            {
                Param* toBeMergedInto = this->At( it->first);
                if( toBeMergedInto->IsNode() && it->second->IsNode() )
                {
                    toBeMergedInto->AsNode().Merge( it->second->AsNode() );
                }
                else
                {
                    this->Replace( it->first, *it->second );
                }
            }
            else
            {
                this->Add( it->first, *it->second );
            }
        }
        return;
    }

    void ParamNode::Erase( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it != fContents.end() )
        {
            delete it->second;
            fContents.erase( it );
        }
        return;
    }

    Param* ParamNode::Remove( const std::string& aName )
    {
        iterator it = fContents.find( aName );
        if( it != fContents.end() )
        {
            Param* removed = it->second;
            fContents.erase( it );
            return removed;
        }
        return NULL;
    }

    ParamNode::iterator ParamNode::Begin()
    {
        return fContents.begin();
    }

    ParamNode::const_iterator ParamNode::Begin() const
    {
        return fContents.begin();
    }

    ParamNode::iterator ParamNode::End()
    {
        return fContents.end();
    }

    ParamNode::const_iterator ParamNode::End() const
    {
        return fContents.end();
    }

    std::string ParamNode::ToString() const
    {
        stringstream out;
        string indentation;
        for ( unsigned i=0; i<Param::sIndentLevel; ++i )
            indentation += "    ";
        out << '\n' << indentation << "{\n";
        Param::sIndentLevel++;
        for( const_iterator it = Begin(); it != End(); ++it )
        {
            out << indentation << "    " << it->first << " : " << *(it->second) << '\n';
        }
        Param::sIndentLevel--;
        out << indentation << "}";
        return out.str();
    }




    std::ostream& operator<<(std::ostream& out, const Param& aValue)
    {
        return out << aValue.ToString();
    }


    std::ostream& operator<<(std::ostream& out, const ParamValue& aValue)
    {
        return out << aValue.ToString();
    }


    std::ostream& operator<<(std::ostream& out, const ParamArray& aValue)
    {
        return out << aValue.ToString();
    }


    std::ostream& operator<<(std::ostream& out, const ParamNode& aValue)
    {
        return out << aValue.ToString();
    }





    ParamInputJSON::ParamInputJSON()
    {
    }
    ParamInputJSON::~ParamInputJSON()
    {
    }

    ParamNode* ParamInputJSON::ReadFile( const std::string& aFilename )
    {
        FILE* t_config_file = fopen( aFilename.c_str(), "r" );
        if( t_config_file == NULL )
        {
            LMCERROR( plog, "file <" << aFilename << "> did not open" );
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
                LMCERROR( plog, "error parsing config file :\n" <<
                        '\t' << t_config_doc.GetParseError() << '\n' <<
                        "\tThe error was reported at line " << newlineCount << ", character " << errorPos - lastNewlinePos );
            }
            else
            {
                LMCERROR( plog, "error parsing config file :\n" <<
                        '\t' << t_config_doc.GetParseError() <<
                        "\tend of file reached before error location was found" );
            }
            fclose( t_config_file );
            return NULL;
        }
        fclose( t_config_file );

        return ParamInputJSON::ReadDocument( t_config_doc );
    }

    ParamNode* ParamInputJSON::ReadString( const std::string& aJSONString )
    {
        rapidjson::Document t_config_doc;
        if( t_config_doc.Parse<0>( aJSONString.c_str() ).HasParseError() )
        {
            LMCERROR( plog, "error parsing string:\n" << t_config_doc.GetParseError() );
            return NULL;
        }
        return ParamInputJSON::ReadDocument( t_config_doc );
    }

    ParamNode* ParamInputJSON::ReadDocument( const rapidjson::Document& aDoc )
    {
        ParamNode* t_config = new ParamNode();
        for( rapidjson::Value::ConstMemberIterator jsonIt = aDoc.MemberBegin();
                jsonIt != aDoc.MemberEnd();
                ++jsonIt)
        {
            t_config->Replace( jsonIt->name.GetString(), ParamInputJSON::ReadValue( jsonIt->value ) );
        }
        return t_config;
    }

    Param* ParamInputJSON::ReadValue( const rapidjson::Value& aValue )
    {
        if( aValue.IsNull() )
        {
            return new Param();
        }
        if( aValue.IsObject() )
        {
            ParamNode* t_config_object = new ParamNode();
            for( rapidjson::Value::ConstMemberIterator jsonIt = aValue.MemberBegin();
                    jsonIt != aValue.MemberEnd();
                    ++jsonIt)
            {
                t_config_object->Replace( jsonIt->name.GetString(), ParamInputJSON::ReadValue( jsonIt->value ) );
            }
            return t_config_object;
        }
        if( aValue.IsArray() )
        {
            ParamArray* t_config_array = new ParamArray();
            for( rapidjson::Value::ConstValueIterator jsonIt = aValue.Begin();
                    jsonIt != aValue.End();
                    ++jsonIt)
            {
                t_config_array->PushBack( ParamInputJSON::ReadValue( *jsonIt ) );
            }
            return t_config_array;
        }
        if( aValue.IsString() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetString();
            return t_config_value;
        }
        if( aValue.IsBool() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetBool();
            return t_config_value;
        }
        if( aValue.IsInt() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetInt();
            return t_config_value;
        }
        if( aValue.IsUint() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetUint();
            return t_config_value;
        }
        if( aValue.IsInt64() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetInt64();
            return t_config_value;
        }
        if( aValue.IsUint64() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetUint64();
            return t_config_value;
        }
        if( aValue.IsDouble() )
        {
            ParamValue* t_config_value = new ParamValue();
            (*t_config_value) << aValue.GetDouble();
            return t_config_value;
        }
        LMCWARN( plog, "(config_reader_json) unknown type; returning null value" );
        return new Param();
    }



    ParamOutputJSON::ParamOutputJSON()
    {}

    ParamOutputJSON::~ParamOutputJSON()
    {}

    bool ParamOutputJSON::WriteFile( const Param& toWrite, const std::string& aFilename, JSONWritingStyle aStyle )
    {
        if( aFilename.empty() )
        {
            LMCERROR( plog, "Filename cannot be an empty string" );
            return false;
        }

        FILE* file = fopen( aFilename.c_str(), "w" );
        if( file == NULL )
        {
            LMCERROR( plog, "Unable to open file: " << aFilename );
            return false;
        }

        rapidjson::FileStream fileStream( file );

        RJWriter* writer = NULL;
        if( aStyle == kCompact )
        {
            writer = new RJWriter( fileStream );
        }
        else
        {
            writer = new RJPrettyWriter( fileStream );
        }

        if (! ParamOutputJSON::WriteParam( toWrite, writer ) )
        {
            LMCERROR( plog, "Error while writing file" );
            delete writer;
            return false;
        }

        delete writer;

        return true;
    }
    bool ParamOutputJSON::WriteParam( const Param& toWrite, RJWriter* writer )
    {
        writer->Null();
        LMCWARN( plog, "writing null" );
        return true;
    }
    bool ParamOutputJSON::WriteParam( const ParamValue& toWrite, RJWriter* writer )
    {
        writer->String(toWrite.ToString().c_str());
        return true;
    }
    bool ParamOutputJSON::WriteParam( const ParamArray& toWrite, RJWriter* writer )
    {
        writer->StartArray();
        for( ParamArray::const_iterator it = toWrite.Begin(); it != toWrite.End(); ++it )
        {
            if( ! ParamOutputJSON::WriteParam( *(*it), writer ) )
            {
                LMCERROR( plog, "Error while writing parameter array" );
                return false;
            }
        }
        writer->EndArray();
        return true;
    }
    bool ParamOutputJSON::WriteParam( const ParamNode& toWrite, RJWriter* writer )
    {
        writer->StartObject();
        for( ParamNode::const_iterator it = toWrite.Begin(); it != toWrite.End(); ++it )
        {
            writer->String( it->first.c_str() );
            if( ! ParamOutputJSON::WriteParam( *(it->second), writer ) )
            {
                LMCERROR( plog, "Error while writing parameter node" );
                return false;
            }
        }
        writer->EndObject();
        return true;
    }


} /* namespace locust */
