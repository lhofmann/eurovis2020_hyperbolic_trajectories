#ifndef SIMPLE_SYSTEM_ADAPTER_H_
#define SIMPLE_SYSTEM_ADAPTER_H_

namespace simple {

template <typename Field>
struct SystemAdapter {
  typedef Eigen::Vector2d argument_type;
  typedef typename Field::value_type value_type;

  Field field;

  explicit SystemAdapter(const Field& field_) : field(field_) {}

  int dimension() const { return 2; }
  value_type zero() const { return field.zero(); }

  bool in_domain(const argument_type &x, double t) const {
    Eigen::Vector3d x_hat;
    x_hat.template head<2>() = x;
    x_hat[2] = t;
    return field.in_domain(x_hat);
  }

  bool operator()(const argument_type &x, double t, Eigen::Ref<value_type> value) {
    Eigen::Vector3d x_hat;
    x_hat.template head<2>() = x;
    x_hat[2] = t;
    return field(x_hat, value);
  }
};

template <typename Field>
SystemAdapter<Field> system_adapter(const Field& field) {
  return SystemAdapter<Field>(field);
}

}

#endif //SIMPLE_SYSTEM_ADAPTER_H_
