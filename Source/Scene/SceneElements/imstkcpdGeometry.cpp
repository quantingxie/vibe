#include "imstkcpdGeometry.h"

namespace cpd
{

  void Geometry::print() const
  {
    std::cout << "Position: " << m_position << ", Orientation: " << m_orientation << std::endl;
  }

  void cpd::Plane::setLengthRef(const Vec3d& p_lengthRef)
  {
    m_lengthRef = p_lengthRef;
    Vec3d n = m_orientation.cross(m_lengthRef);

    if (n.norm() == 0)
    {
      unsigned i = 0;
      std::cout << "LengthRef aligns with geometry orientation!" << std::endl;

      Vec3d axis[3];
      axis[0] = CPD_UP_VECTOR;  axis[1] = CPD_RIGHT_VECTOR;  axis[2] = CPD_FORWARD_VECTOR;
      while (n.norm() == 0)
      {
        n = m_orientation.cross(axis[i]);
        i++;
      }

      std::cout << "New lengthRef is set to be " << axis[i] << ", provide new valid lengthRef otherwise." << std::endl;
    }

    m_lengthDirection = n.normalized();
    m_widthDirection = m_orientation.cross(m_lengthDirection);
  }

  void cpd::Plane::update()
  {
    m_keyVertices[0] = m_position - m_lengthDirection - m_widthDirection;
    m_keyVertices[1] = m_position + m_lengthDirection - m_widthDirection;
    m_keyVertices[2] = m_position + m_lengthDirection + m_widthDirection;
    m_keyVertices[3] = m_position - m_lengthDirection + m_widthDirection;
  }

  void Plane::print() const
  {
    std::cout << "Type: Plane. " << std::endl;
    Geometry::print();
    std::cout << "Length: " << m_halfLength * 2 << "Width: " << m_halfWidth * 2 << std::endl;
  }

  const double Plane::getVolume() const
  {
    return 0.0;
  }



  void Sphere::update()
  {
    m_keyVertices[0] = m_position;
  }

  void Sphere::print() const
  {
    std::cout << "Type: Sphere. " << std::endl;
    Geometry::print();
    std::cout << "Radius: " << m_radius << std::endl;
  }

  const double Sphere::getVolume() const
  {
    return 4.0 / 3.0* CPD_PI * m_radius*m_radius*m_radius;
  }


  void Cylinder::update()
  {
    m_keyVertices[0] = m_position;
    m_keyVertices[1] = m_position + m_orientation * m_length;
  }

  void Cylinder::print() const
  {
    std::cout << "Type: Cylinder. " << std::endl;
    Geometry::print();
    std::cout << "Radius: " << m_radius << ", Length: " << m_length << std::endl;
  }

  const double Cylinder::getVolume() const
  {
    return CPD_PI * m_radius*m_radius*m_length;
  }

}