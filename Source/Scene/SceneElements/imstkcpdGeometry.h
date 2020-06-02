#ifndef CPDGEOMETRY_H
#define CPDGEOMETRY_H

#include <memory>

#include "imstkcpdMath.h"

namespace cpd
{
  class Geometry
  {
  public:
    enum class Type
    {
      Plane,
      Sphere,
      Cylinder,
      Cube,
      Undefined
    };

  public:
    Geometry() = delete;
    Geometry(Type p_type) : m_type(p_type) {}
    Geometry(Type p_type, Vec3d p_position, Vec3d p_orientation) : m_type(p_type), m_position(p_position), m_orientation(p_orientation) {}

    const Type getType() const { return m_type; }
    const StdVectorOfVec3d getKeyVertices() const { return m_keyVertices; }

    virtual void update() = 0;

    virtual void print() const;
    virtual const double getVolume() const = 0;

    void setPosition(const Vec3d& p_position) { m_position = p_position; }
    void setOrientation(const Vec3d& p_orientation) { m_orientation = p_orientation; }

    const Vec3d& getPosistion() const { return m_position; }
    const Vec3d& getOrientation() const { return m_orientation; }

  protected:
    Type m_type = Type::Undefined;
    StdVectorOfVec3d m_keyVertices;
    unsigned m_keyVsize;

    Vec3d m_position = CPD_WORLD_ORIGIN;
    Vec3d m_orientation = CPD_UP_VECTOR;
  };

  using GeometryPtr = std::shared_ptr<Geometry>;


  class Plane :public Geometry
  {
  public:
    Plane() : Geometry(Type::Plane) { m_keyVsize = 4; m_keyVertices.resize(m_keyVsize); }
    Plane(double p_length, double p_width) : Geometry(Type::Plane), m_halfLength(p_length / 2.0), m_halfWidth(p_width / 2.0) { m_keyVsize = 4; m_keyVertices.resize(m_keyVsize); }

    void setLengthRef(const Vec3d& p_lengthRef);
    void setLength(double p_length) { m_halfLength = 0.5*p_length; }
    void setWidth(double p_width) { m_halfWidth = 0.5*p_width; }

    const Vec3d& getLengthRef() const { return m_lengthRef; }
    const double getLength() const { return 2.0*m_halfLength; }
    const double getWidth() const { return 2.0*m_halfWidth; }

    void update() override;
    void print() const override;
    const double getVolume() const override;

  private:
    double m_halfLength = 2.0;
    double m_halfWidth = 1.5;
    Vec3d m_lengthRef = CPD_RIGHT_VECTOR;
    Vec3d m_lengthDirection;
    Vec3d m_widthDirection;
  };
  using PlanePtr = std::shared_ptr<Plane>;


  class Sphere :public Geometry
  {
  public:
    Sphere() : Geometry(Type::Sphere) { m_keyVsize = 1; m_keyVertices.resize(m_keyVsize); }
    Sphere(double p_radius) : Geometry(Type::Sphere), m_radius(p_radius) { m_keyVsize = 1; m_keyVertices.resize(m_keyVsize); }

    void update() override;
    void print() const override;
    const double getVolume() const override;

    void setRadius(double p_radius) { m_radius = p_radius; }
    const double getRadius() const { return m_radius; }

  private:
    double m_radius = 1.0;
  };
  using SpherePtr = std::shared_ptr<Sphere>;


  class Cylinder :public Geometry
  {
  public:
    Cylinder() :Geometry(Type::Cylinder) { m_keyVsize = 2; m_keyVertices.resize(m_keyVsize); }
    Cylinder(double p_radius, double p_length) : Geometry(Type::Cylinder), m_radius(p_radius), m_length(p_length) { 
      m_keyVsize = 2; m_keyVertices.resize(m_keyVsize); 
    }
    Cylinder(double p_radius, double p_length, Vec3d p_position, Vec3d p_orientation): Geometry(Type::Cylinder, p_position, p_orientation),
      m_radius(p_radius), m_length(p_length) {
      m_keyVsize = 2; m_keyVertices.resize(m_keyVsize);
    }

    void update() override;
    void print() const override;
    const double getVolume() const override;

    void setRadius(double p_radius) { m_radius = p_radius; }
    void setLength(double p_length) { m_radius = p_length; }

    const double getRadius() const { return m_radius; }
    const double getLength() const { return m_length; }

  private:
    double m_radius = 1.0;
    double m_length = 5.0;
  };
  using CylinderPtr = std::shared_ptr<Cylinder>;

}


#endif // !CPDGEOMETRY_H
