// Code generated by Stan version 2.14

#ifndef VIFCOPULA_SERVICE_WRITE_THETA_HPP
#define VIFCOPULA_SERVICE_WRITE_THETA_HPP

namespace vifcopula{

void write_theta(int id, int copula_type,
                stan::io::reader<double>& in__,
                std::vector<double>& vars__){
    double theta;
    double theta2;

    // read-transform, write parameters
    switch ( copula_type )
    {
    case 0:
        // Independent copula
        break;
    case 1:
        // Gaussian copula
        // if (id == 0){
        //     theta = in__.scalar_lub_constrain(0,1);
        // } else {
        //     theta = in__.scalar_lub_constrain(-1,1);
        // }
        theta = in__.scalar_lub_constrain(-1,1);
        vars__.push_back(theta);
        break;
    case 2:
        // Student copula
        theta = in__.scalar_lub_constrain(0,1);
        theta2 = in__.scalar_lub_constrain(2,15);
        vars__.push_back(theta);
        vars__.push_back(theta2);
        break;
    case 3:
        // Clayon copula
        theta = in__.scalar_lub_constrain(0.001,30);
        vars__.push_back(theta);
        break;
    case 4:
        // Gumbel copula
        theta = in__.scalar_lub_constrain(1,15);
        vars__.push_back(theta);
        break;
    case 5:
        // Frank copula
        theta = in__.scalar_lub_constrain(0,100);
        vars__.push_back(theta);
        break;
    case 6:
        // Joe copula
        theta = in__.scalar_lub_constrain(1,30);
        vars__.push_back(theta);
        break;
    default:
        // Code to execute if <variable> does not equal the value following any of the cases
        // Send a break message.
        break;
    }


}

} // namespace
#endif // VIFCOPULA_BICOPULA_HPP