/*************************************************************************/
/*  register_types.cpp                                                   */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2019 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2019 Godot Engine contributors (cf. AUTHORS.md)    */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "register_types.h"

#include "godot_module_files/bullet_physics_server.h"
#include "core/class_db.h"
#include "core/project_settings.h"

#include <iostream>

/**
	@author AndreaCatania
*/

#ifndef _3D_DISABLED
PhysicsServer *_createCustomPhysicsCallback() {
	return memnew(CustomPhysicsServer);
}
#endif

void register_customphysics_types() {
#ifndef _3D_DISABLED
	PhysicsServerManager::register_server("CustomPhysics", &_createCustomPhysicsCallback);

	GLOBAL_DEF("physics/3d/active_soft_world", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/3d/active_soft_world",
			PropertyInfo(Variant::BOOL, "physics/3d/active_soft_world"));

	GLOBAL_DEF("physics/customphysics/constraint_iterations", 30);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/constraint_iterations",
			PropertyInfo(Variant::INT, "physics/customphysics/constraint_iterations"));

	GLOBAL_DEF("physics/customphysics/hinge_iterations", 10);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/hinge_iterations",
			PropertyInfo(Variant::INT, "physics/customphysics/hinge_iterations"));

	GLOBAL_DEF("physics/customphysics/contact_gamma", 0.1f);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/contact_gamma",
			PropertyInfo(Variant::REAL, "physics/customphysics/contact_gamma"));

	GLOBAL_DEF("physics/customphysics/ball_joint_gamma", 0.05f);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/ball_joint_gamma",
			PropertyInfo(Variant::REAL, "physics/customphysics/ball_joint_gamma"));

	GLOBAL_DEF("physics/customphysics/hinge_joint_gamma", 0.1f);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/hinge_joint_gamma",
			PropertyInfo(Variant::REAL, "physics/customphysics/hinge_joint_gamma"));

	GLOBAL_DEF("physics/customphysics/mu", 0.3f);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/mu",
			PropertyInfo(Variant::REAL, "physics/customphysics/mu"));

	GLOBAL_DEF("physics/customphysics/friction_constraints", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/friction_constraints",
			PropertyInfo(Variant::BOOL, "physics/customphysics/friction_constraints"));
	
	GLOBAL_DEF("physics/customphysics/contact_constraints", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/contact_constraints",
			PropertyInfo(Variant::BOOL, "physics/customphysics/contact_constraints"));
	
	GLOBAL_DEF("physics/customphysics/ball_joints_constraints", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/ball_joints_constraints",
			PropertyInfo(Variant::BOOL, "physics/customphysics/ball_joints_constraints"));

	GLOBAL_DEF("physics/customphysics/hinge_joints_constraints", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/hinge_joints_constraints",
			PropertyInfo(Variant::BOOL, "physics/customphysics/hinge_joints_constraints"));
			
	GLOBAL_DEF("physics/customphysics/warm_starting", true);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/warm_starting",
			PropertyInfo(Variant::BOOL, "physics/customphysics/warm_starting"));

	GLOBAL_DEF("physics/customphysics/warm_starting_factor", 1.f);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/warm_starting_factor",
			PropertyInfo(Variant::REAL, "physics/customphysics/warm_starting_factor"));

	GLOBAL_DEF("physics/customphysics/hinge_with_2x2", false);
	ProjectSettings::get_singleton()->set_custom_property_info("physics/customphysics/hinge_with_2x2",
			PropertyInfo(Variant::BOOL, "physics/customphysics/hinge_with_2x2"));
#endif
}

void unregister_customphysics_types() {
}
