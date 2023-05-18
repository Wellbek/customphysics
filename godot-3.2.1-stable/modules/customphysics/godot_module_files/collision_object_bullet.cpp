/*************************************************************************/
/*  collision_object_bullet.cpp                                          */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2020 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2020 Godot Engine contributors (cf. AUTHORS.md).   */
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

#include "collision_object_bullet.h"

#include "area_bullet.h"
#include "bullet_physics_server.h"
#include "bullet_types_converter.h"
#include "bullet_utilities.h"
#include "shape_bullet.h"
#include "space_bullet.h"

#include <btBulletCollisionCommon.h>

/**
	@author AndreaCatania
*/

// We enable dynamic AABB tree so that we can actually perform a broadphase on bodies with compound collision shapes.
// This is crucial for the performance of kinematic bodies and for bodies with transforming shapes.
#define enableDynamicAabbTree true

CollisionObjectCustom::ShapeWrapper::~ShapeWrapper() {}

void CollisionObjectCustom::ShapeWrapper::set_transform(const Transform &p_transform) {
	G_TO_B_CUSTOM(p_transform.get_basis().get_scale_abs(), scale);
	G_TO_B_CUSTOM(p_transform, transform);
	UNSCALE_BT_BASIS_CUSTOM(transform);
}

void CollisionObjectCustom::ShapeWrapper::set_transform(const btTransform &p_transform) {
	transform = p_transform;
}

btTransform CollisionObjectCustom::ShapeWrapper::get_adjusted_transform() const {
	if (shape->get_type() == PhysicsServer::SHAPE_HEIGHTMAP) {
		const HeightMapShapeCustom *hm_shape = (const HeightMapShapeCustom *)shape; // should be safe to cast now
		btTransform adjusted_transform;

		// Bullet centers our heightmap:
		// https://github.com/bulletphysics/bullet3/blob/master/src/BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h#L33
		// This is really counter intuitive so we're adjusting for it

		adjusted_transform.setIdentity();
		adjusted_transform.setOrigin(btVector3(0.0, hm_shape->min_height + ((hm_shape->max_height - hm_shape->min_height) * 0.5), 0.0));
		adjusted_transform *= transform;

		return adjusted_transform;
	} else {
		return transform;
	}
}

void CollisionObjectCustom::ShapeWrapper::claim_bt_shape(const btVector3 &body_scale) {
	if (!bt_shape) {
		if (active)
			bt_shape = shape->create_bt_shape(scale * body_scale);
		else
			bt_shape = ShapeCustom::create_shape_empty();
	}
}

CollisionObjectCustom::CollisionObjectCustom(Type p_type) :
		RIDCustom(),
		type(p_type),
		instance_id(0),
		collisionLayer(0),
		collisionMask(0),
		collisionsEnabled(true),
		m_isStatic(false),
		ray_pickable(false),
		bt_collision_object(NULL),
		body_scale(1., 1., 1.),
		force_shape_reset(false),
		space(NULL),
		isTransformChanged(false) {}

CollisionObjectCustom::~CollisionObjectCustom() {
	// Remove all overlapping, notify is not required since godot take care of it
	for (int i = areasOverlapped.size() - 1; 0 <= i; --i) {
		areasOverlapped[i]->remove_overlap(this, /*Notify*/ false);
	}

	destroyCustomCollisionObject();
}

bool equal_custom(real_t first, real_t second) {
	return Math::abs(first - second) <= 0.001f;
}

void CollisionObjectCustom::set_body_scale(const Vector3 &p_new_scale) {
	if (!equal_custom(p_new_scale[0], body_scale[0]) || !equal_custom(p_new_scale[1], body_scale[1]) || !equal_custom(p_new_scale[2], body_scale[2])) {
		body_scale = p_new_scale;
		body_scale_changed();
	}
}

btVector3 CollisionObjectCustom::get_bt_body_scale() const {
	btVector3 s;
	G_TO_B_CUSTOM(body_scale, s);
	return s;
}

void CollisionObjectCustom::body_scale_changed() {
	force_shape_reset = true;
}

void CollisionObjectCustom::destroyCustomCollisionObject() {
	bulletdelete(bt_collision_object);
}

void CollisionObjectCustom::setupCustomCollisionObject(btCollisionObject *p_collisionObject) {
	bt_collision_object = p_collisionObject;
	bt_collision_object->setUserPointer(this);
	bt_collision_object->setUserIndex(type);
	// Force the enabling of collision and avoid problems
	set_collision_enabled(collisionsEnabled);
	p_collisionObject->setCollisionFlags(p_collisionObject->getCollisionFlags() | btCollisionObject::CF_CUSTOM_MATERIAL_CALLBACK);
}

void CollisionObjectCustom::add_collision_exception(const CollisionObjectCustom *p_ignoreCollisionObject) {
	exceptions.insert(p_ignoreCollisionObject->get_self());
	if (!bt_collision_object)
		return;
	bt_collision_object->setIgnoreCollisionCheck(p_ignoreCollisionObject->bt_collision_object, true);
	if (space)
		space->get_broadphase()->getOverlappingPairCache()->cleanProxyFromPairs(bt_collision_object->getBroadphaseHandle(), space->get_dispatcher());
}

void CollisionObjectCustom::remove_collision_exception(const CollisionObjectCustom *p_ignoreCollisionObject) {
	exceptions.erase(p_ignoreCollisionObject->get_self());
	bt_collision_object->setIgnoreCollisionCheck(p_ignoreCollisionObject->bt_collision_object, false);
	if (space)
		space->get_broadphase()->getOverlappingPairCache()->cleanProxyFromPairs(bt_collision_object->getBroadphaseHandle(), space->get_dispatcher());
}

bool CollisionObjectCustom::has_collision_exception(const CollisionObjectCustom *p_otherCollisionObject) const {
	return !bt_collision_object->checkCollideWith(p_otherCollisionObject->bt_collision_object);
}

void CollisionObjectCustom::set_collision_enabled(bool p_enabled) {
	collisionsEnabled = p_enabled;
	if (collisionsEnabled) {
		bt_collision_object->setCollisionFlags(bt_collision_object->getCollisionFlags() & (~btCollisionObject::CF_NO_CONTACT_RESPONSE));
	} else {
		bt_collision_object->setCollisionFlags(bt_collision_object->getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
	}
}

bool CollisionObjectCustom::is_collisions_response_enabled() {
	return collisionsEnabled;
}

void CollisionObjectCustom::notify_new_overlap(AreaCustom *p_area) {
	areasOverlapped.push_back(p_area);
}

void CollisionObjectCustom::on_exit_area(AreaCustom *p_area) {
	areasOverlapped.erase(p_area);
}

void CollisionObjectCustom::set_godot_object_flags(int flags) {
	bt_collision_object->setUserIndex2(flags);
}

int CollisionObjectCustom::get_godot_object_flags() const {
	return bt_collision_object->getUserIndex2();
}

void CollisionObjectCustom::set_transform(const Transform &p_global_transform) {

	set_body_scale(p_global_transform.basis.get_scale_abs());

	btTransform bt_transform;
	G_TO_B_CUSTOM(p_global_transform, bt_transform);
	UNSCALE_BT_BASIS_CUSTOM(bt_transform);

	set_transform__bullet(bt_transform);
}

Transform CollisionObjectCustom::get_transform() const {
	Transform t;
	B_TO_G_CUSTOM(get_transform__bullet(), t);
	t.basis.scale(body_scale);
	return t;
}

void CollisionObjectCustom::set_transform__bullet(const btTransform &p_global_transform) {
	bt_collision_object->setWorldTransform(p_global_transform);
	notify_transform_changed();
}

const btTransform &CollisionObjectCustom::get_transform__bullet() const {
	return bt_collision_object->getWorldTransform();
}

void CollisionObjectCustom::notify_transform_changed() {
	isTransformChanged = true;
}

RigidCollisionObjectCustom::RigidCollisionObjectCustom(Type p_type) :
		CollisionObjectCustom(p_type),
		mainShape(NULL) {
}

RigidCollisionObjectCustom::~RigidCollisionObjectCustom() {
	remove_all_shapes(true, true);
	if (mainShape && mainShape->isCompound()) {
		bulletdelete(mainShape);
	}
}

void RigidCollisionObjectCustom::add_shape(ShapeCustom *p_shape, const Transform &p_transform, bool p_disabled) {
	shapes.push_back(ShapeWrapper(p_shape, p_transform, !p_disabled));
	p_shape->add_owner(this);
	reload_shapes();
}

void RigidCollisionObjectCustom::set_shape(int p_index, ShapeCustom *p_shape) {
	ShapeWrapper &shp = shapes.write[p_index];
	shp.shape->remove_owner(this);
	p_shape->add_owner(this);
	shp.shape = p_shape;
	reload_shapes();
}

int RigidCollisionObjectCustom::get_shape_count() const {
	return shapes.size();
}

ShapeCustom *RigidCollisionObjectCustom::get_shape(int p_index) const {
	return shapes[p_index].shape;
}

btCollisionShape *RigidCollisionObjectCustom::get_bt_shape(int p_index) const {
	return shapes[p_index].bt_shape;
}

int RigidCollisionObjectCustom::find_shape(ShapeCustom *p_shape) const {
	const int size = shapes.size();
	for (int i = 0; i < size; ++i) {
		if (shapes[i].shape == p_shape)
			return i;
	}
	return -1;
}

void RigidCollisionObjectCustom::remove_shape_full(ShapeCustom *p_shape) {
	// Remove the shape, all the times it appears
	// Reverse order required for delete.
	for (int i = shapes.size() - 1; 0 <= i; --i) {
		if (p_shape == shapes[i].shape) {
			internal_shape_destroy(i);
			shapes.remove(i);
		}
	}
	reload_shapes();
}

void RigidCollisionObjectCustom::remove_shape_full(int p_index) {
	ERR_FAIL_INDEX(p_index, get_shape_count());
	internal_shape_destroy(p_index);
	shapes.remove(p_index);
	reload_shapes();
}

void RigidCollisionObjectCustom::remove_all_shapes(bool p_permanentlyFromThisBody, bool p_force_not_reload) {
	// Reverse order required for delete.
	for (int i = shapes.size() - 1; 0 <= i; --i) {
		internal_shape_destroy(i, p_permanentlyFromThisBody);
	}
	shapes.clear();
	if (!p_force_not_reload)
		reload_shapes();
}

void RigidCollisionObjectCustom::set_shape_transform(int p_index, const Transform &p_transform) {
	ERR_FAIL_INDEX(p_index, get_shape_count());

	shapes.write[p_index].set_transform(p_transform);
	shape_changed(p_index);
}

const btTransform &RigidCollisionObjectCustom::get_bt_shape_transform(int p_index) const {
	return shapes[p_index].transform;
}

Transform RigidCollisionObjectCustom::get_shape_transform(int p_index) const {
	Transform trs;
	B_TO_G_CUSTOM(shapes[p_index].transform, trs);
	return trs;
}

void RigidCollisionObjectCustom::set_shape_disabled(int p_index, bool p_disabled) {
	if (shapes[p_index].active != p_disabled)
		return;
	shapes.write[p_index].active = !p_disabled;
	shape_changed(p_index);
}

bool RigidCollisionObjectCustom::is_shape_disabled(int p_index) {
	return !shapes[p_index].active;
}

void RigidCollisionObjectCustom::shape_changed(int p_shape_index) {
	ShapeWrapper &shp = shapes.write[p_shape_index];
	if (shp.bt_shape == mainShape) {
		mainShape = NULL;
	}
	bulletdelete(shp.bt_shape);
	reload_shapes();
}

void RigidCollisionObjectCustom::reload_shapes() {

	if (mainShape && mainShape->isCompound()) {
		// Destroy compound
		bulletdelete(mainShape);
	}

	mainShape = NULL;

	ShapeWrapper *shpWrapper;
	const int shape_count = shapes.size();

	// Reset shape if required
	if (force_shape_reset) {
		for (int i(0); i < shape_count; ++i) {
			shpWrapper = &shapes.write[i];
			bulletdelete(shpWrapper->bt_shape);
		}
		force_shape_reset = false;
	}

	const btVector3 body_scale(get_bt_body_scale());

	// Try to optimize by not using compound
	if (1 == shape_count) {
		shpWrapper = &shapes.write[0];
		btTransform transform = shpWrapper->get_adjusted_transform();
		if (transform.getOrigin().isZero() && transform.getBasis() == transform.getBasis().getIdentity()) {
			shpWrapper->claim_bt_shape(body_scale);
			mainShape = shpWrapper->bt_shape;
			main_shape_changed();
			return;
		}
	}

	// Optimization not possible use a compound shape
	btCompoundShape *compoundShape = bulletnew(btCompoundShape(enableDynamicAabbTree, shape_count));

	for (int i(0); i < shape_count; ++i) {
		shpWrapper = &shapes.write[i];
		shpWrapper->claim_bt_shape(body_scale);
		btTransform scaled_shape_transform(shpWrapper->get_adjusted_transform());
		scaled_shape_transform.getOrigin() *= body_scale;
		compoundShape->addChildShape(scaled_shape_transform, shpWrapper->bt_shape);
	}

	compoundShape->recalculateLocalAabb();
	mainShape = compoundShape;
	main_shape_changed();
}

void RigidCollisionObjectCustom::body_scale_changed() {
	CollisionObjectCustom::body_scale_changed();
	reload_shapes();
}

void RigidCollisionObjectCustom::internal_shape_destroy(int p_index, bool p_permanentlyFromThisBody) {
	ShapeWrapper &shp = shapes.write[p_index];
	shp.shape->remove_owner(this, p_permanentlyFromThisBody);
	if (shp.bt_shape == mainShape) {
		mainShape = NULL;
	}
	bulletdelete(shp.bt_shape);
}