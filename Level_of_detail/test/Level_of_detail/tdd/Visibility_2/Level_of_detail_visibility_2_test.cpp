// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

// New CGAL includes.
#include <CGAL/Visibility_2/Level_of_detail_visibility_2.h>
#include <CGAL/Mylog/Mylog.h>

using namespace testing;

class LOD_VisibilityTest: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point_2   = Traits::Point_2;
	using Point_3   = Traits::Point_3;
	using Container = CGAL::Point_set_3<Point_3>;
	using Iterator  = Container::iterator;

	using Label     = int; 
	using Label_map = Container:: template Property_map<Label>; 

	using VB = CGAL::Triangulation_vertex_base_2<Traits>;
	using FB = CGAL::Constrained_triangulation_face_base_2<Traits>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using LodVisibility = CGAL::LOD::Level_of_detail_visibility_from_classification_2<Traits, Container, CDT>;

	using Vertex_handle     = CDT::Vertex_handle;
	using Visibility_result = std::map<int, LodVisibility::Visibility_label>;

	using Log = CGAL::LOD::Mylog;

	using Segment  = Traits::Segment_2;
	using Segments = std::vector<Segment>;

	LodVisibility lodVisibility;
	CDT cdt; Container input;

	LOD_VisibilityTest() {
		create_data();
	}

	void create_data() {

		cdt.clear();
		input.clear();		

		set_basic_input(cdt, input);
	}

	void set_labels_property(Container &input, Label_map &labels) {
		
		auto success = false;
		boost::tie(labels, success)  = input. template add_property_map<Label>("label", -1);
		assert(success);
	}

	void set_basic_input(CDT &cdt, Container &input) {

		Label_map labels;
		set_labels_property(input, labels);

		const Label ground     = 0;
		const Label facade     = 1;
		const Label roof       = 2;
		const Label vegetation = 3; 

		Iterator 
		it = input.insert(Point_3(0.25, 0.25, 0)); labels[*it] = roof;
		it = input.insert(Point_3(0.70, 0.70, 0)); labels[*it] = roof;
		it = input.insert(Point_3(1.30, 0.30, 0)); labels[*it] = ground;
		it = input.insert(Point_3(1.60, 0.60, 0)); labels[*it] = vegetation;
		it = input.insert(Point_3(1.00, 0.80, 0)); labels[*it] = facade;

		Vertex_handle va = cdt.insert(Point_2(0, 0));
		Vertex_handle vb = cdt.insert(Point_2(1, 0));
		Vertex_handle vc = cdt.insert(Point_2(0, 1));
		Vertex_handle vd = cdt.insert(Point_2(1, 1));
		Vertex_handle ve = cdt.insert(Point_2(2, 0));
		Vertex_handle vf = cdt.insert(Point_2(2, 1));

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, va);

		cdt.insert_constraint(vb, vd);
		cdt.insert_constraint(vd, vc);
		cdt.insert_constraint(vc, vb);

		cdt.insert_constraint(vb, ve);
		cdt.insert_constraint(ve, vd);
		cdt.insert_constraint(vd, vb);

		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vd);
		cdt.insert_constraint(vd, ve);		
	}
};

TEST_F(LOD_VisibilityTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_VisibilityTest, SavesVisibility) {

	Visibility_result visibility; Segments segments;
	lodVisibility.compute(cdt, input, visibility);

	Log log;
	log.save_visibility_eps(cdt, visibility, input, segments);
}

TEST_F(LOD_VisibilityTest, VerifiesLabels) {

	Visibility_result visibility;
	lodVisibility.compute(cdt, input, visibility);

	ASSERT_THAT(visibility.size(), Eq(static_cast<size_t>(4)));

	ASSERT_THAT(visibility[0], Eq(LodVisibility::Visibility_label::IN));
	ASSERT_THAT(visibility[1], Eq(LodVisibility::Visibility_label::IN));

	ASSERT_THAT(visibility[2], Eq(LodVisibility::Visibility_label::OUT));
	ASSERT_THAT(visibility[3], Eq(LodVisibility::Visibility_label::OUT));
}
