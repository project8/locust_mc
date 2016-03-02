/*
 * LMCSingleton.hh
 *
 *  Created on: Jul 31, 2012
 *      Author: nsoblath
 */

#ifndef LMCFACTORY_HH_
#define LMCFACTORY_HH_

#include "LMCSingleton.hh"

#include "Logger.hh"

#include <map>
#include <string>
#include <utility>

namespace locust
{
    LOGGER( lmclog_fact, "Factory" );

    template< class XBaseType >
    class Factory;

    template< class XBaseType >
    class BaseRegistrar
    {
        public:
            BaseRegistrar() {}
            virtual ~BaseRegistrar() {}

        public:
            friend class Factory< XBaseType >;

        protected:
            virtual XBaseType* Create() const = 0;

    };

    template< class XBaseType, class XDerivedType >
    class Registrar : public BaseRegistrar< XBaseType >
    {
        public:
            Registrar(const std::string& a_class_name);
            virtual ~Registrar();

        protected:
            void RegisterClass(const std::string& a_class_name) const;

            XBaseType* Create() const;

    };


    template< class XBaseType >
    class Factory : public Singleton< Factory< XBaseType > >
    {
        public:
            typedef std::map< std::string, const BaseRegistrar< XBaseType >* > FactoryMap;
            typedef typename FactoryMap::value_type FactoryEntry;
            typedef typename FactoryMap::iterator FactoryIt;
            typedef typename FactoryMap::const_iterator FactoryCIt;

        public:
            XBaseType* Create(const std::string& a_class_name);
            XBaseType* Create(const FactoryCIt& iter);

            void RegisterClass(const std::string& a_class_name, const BaseRegistrar< XBaseType >* BaseRegistrar);

        protected:
            FactoryMap* fMap;


        protected:
            friend class Singleton< Factory >;
            friend class Destroyer< Factory >;
            Factory();
            ~Factory();
    };

    template< class XBaseType >
    XBaseType* Factory< XBaseType >::Create(const FactoryCIt& iter)
    {
        return iter->second->Create();
    }

    template< class XBaseType >
    XBaseType* Factory< XBaseType >::Create(const std::string& a_class_name)
    {
        FactoryCIt it = fMap->find(a_class_name);
        if (it == fMap->end())
        {
            ERROR( lmclog_fact, "Did not find Factory for <" << a_class_name << ">." );
            return NULL;
        }

        return it->second->Create();
    }

    template< class XBaseType >
    void Factory< XBaseType >::RegisterClass(const std::string& a_class_name, const BaseRegistrar< XBaseType >* a_Registrar)
    {
        FactoryCIt it = fMap->find(a_class_name);
        if (it != fMap->end())
        {
            ERROR( lmclog_fact, "Already have Factory RegisterClassed for <" << a_class_name << ">." );
            return;
        }
        fMap->insert(std::pair< std::string, const BaseRegistrar< XBaseType >* >(a_class_name, a_Registrar));
        //std::cout << "registered a Factory for class " << a_class_name << ", Factory #" << fMap->size()-1 << std::endl;
    }

    template< class XBaseType >
    Factory< XBaseType >::Factory() :
        fMap(new FactoryMap())
    {}

    template< class XBaseType >
    Factory< XBaseType >::~Factory()
    {
        delete fMap;
    }





    template< class XBaseType, class XDerivedType >
    Registrar< XBaseType, XDerivedType >::Registrar(const std::string& a_class_name) :
            BaseRegistrar< XBaseType >()
    {
        RegisterClass(a_class_name);
    }

    template< class XBaseType, class XDerivedType >
    Registrar< XBaseType, XDerivedType >::~Registrar()
    {}

    template< class XBaseType, class XDerivedType >
    void Registrar< XBaseType, XDerivedType >::RegisterClass(const std::string& a_class_name) const
    {
        Factory< XBaseType >::GetInstance()->RegisterClass(a_class_name, this);
        return;
    }

    template< class XBaseType, class XDerivedType >
    XBaseType* Registrar< XBaseType, XDerivedType >::Create() const
    {
        return dynamic_cast< XBaseType* >(new XDerivedType());
    }

} /* namespace mantis */
#endif /* LMCFACTORY_HH_ */
