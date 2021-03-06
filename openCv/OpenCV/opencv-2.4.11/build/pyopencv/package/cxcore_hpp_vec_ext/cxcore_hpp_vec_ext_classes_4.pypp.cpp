// This file has been generated by Py++.

#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "cxcore_hpp_vec_wrapper.hpp"
#include "opencv_headers.hpp"
#include "cxcore_hpp_vec_ext_classes_4.pypp.hpp"

namespace bp = boost::python;

static inline void vector_Vec3f_resize(::std::vector< cv::Vec<float, 3> > &inst, size_t num) { inst.resize(num); }

static inline void vector_Vec2f_resize(::std::vector< cv::Vec<float, 2> > &inst, size_t num) { inst.resize(num); }

static inline void vector_Vec6d_resize(::std::vector< cv::Vec<double, 6> > &inst, size_t num) { inst.resize(num); }

static inline void vector_Vec4d_resize(::std::vector< cv::Vec<double, 4> > &inst, size_t num) { inst.resize(num); }

static inline void vector_Vec3d_resize(::std::vector< cv::Vec<double, 3> > &inst, size_t num) { inst.resize(num); }

void register_classes_4(){

    { //::std::vector< cv::Vec<float, 3> >
        typedef bp::class_< std::vector< cv::Vec<float, 3> > > vector_Vec3f_exposer_t;
        vector_Vec3f_exposer_t vector_Vec3f_exposer = vector_Vec3f_exposer_t( "vector_Vec3f" );
        bp::scope vector_Vec3f_scope( vector_Vec3f_exposer );
        //WARNING: the next line of code will not compile, because "cv::Vec<float,3>" does not have operator== !
        vector_Vec3f_exposer.def( bp::vector_indexing_suite< ::std::vector< cv::Vec<float, 3> > >() );
        vector_Vec3f_exposer.def("resize", &::vector_Vec3f_resize, ( bp::arg("num") ));
    }

    { //::std::vector< cv::Vec<float, 2> >
        typedef bp::class_< std::vector< cv::Vec<float, 2> > > vector_Vec2f_exposer_t;
        vector_Vec2f_exposer_t vector_Vec2f_exposer = vector_Vec2f_exposer_t( "vector_Vec2f" );
        bp::scope vector_Vec2f_scope( vector_Vec2f_exposer );
        //WARNING: the next line of code will not compile, because "cv::Vec<float,2>" does not have operator== !
        vector_Vec2f_exposer.def( bp::vector_indexing_suite< ::std::vector< cv::Vec<float, 2> > >() );
        vector_Vec2f_exposer.def("resize", &::vector_Vec2f_resize, ( bp::arg("num") ));
    }

    { //::std::vector< cv::Vec<double, 6> >
        typedef bp::class_< std::vector< cv::Vec<double, 6> > > vector_Vec6d_exposer_t;
        vector_Vec6d_exposer_t vector_Vec6d_exposer = vector_Vec6d_exposer_t( "vector_Vec6d" );
        bp::scope vector_Vec6d_scope( vector_Vec6d_exposer );
        //WARNING: the next line of code will not compile, because "cv::Vec<double,6>" does not have operator== !
        vector_Vec6d_exposer.def( bp::vector_indexing_suite< ::std::vector< cv::Vec<double, 6> > >() );
        vector_Vec6d_exposer.def("resize", &::vector_Vec6d_resize, ( bp::arg("num") ));
    }

    { //::std::vector< cv::Vec<double, 4> >
        typedef bp::class_< std::vector< cv::Vec<double, 4> > > vector_Vec4d_exposer_t;
        vector_Vec4d_exposer_t vector_Vec4d_exposer = vector_Vec4d_exposer_t( "vector_Vec4d" );
        bp::scope vector_Vec4d_scope( vector_Vec4d_exposer );
        //WARNING: the next line of code will not compile, because "cv::Vec<double,4>" does not have operator== !
        vector_Vec4d_exposer.def( bp::vector_indexing_suite< ::std::vector< cv::Vec<double, 4> > >() );
        vector_Vec4d_exposer.def("resize", &::vector_Vec4d_resize, ( bp::arg("num") ));
    }

    { //::std::vector< cv::Vec<double, 3> >
        typedef bp::class_< std::vector< cv::Vec<double, 3> > > vector_Vec3d_exposer_t;
        vector_Vec3d_exposer_t vector_Vec3d_exposer = vector_Vec3d_exposer_t( "vector_Vec3d" );
        bp::scope vector_Vec3d_scope( vector_Vec3d_exposer );
        //WARNING: the next line of code will not compile, because "cv::Vec<double,3>" does not have operator== !
        vector_Vec3d_exposer.def( bp::vector_indexing_suite< ::std::vector< cv::Vec<double, 3> > >() );
        vector_Vec3d_exposer.def("resize", &::vector_Vec3d_resize, ( bp::arg("num") ));
    }

}
