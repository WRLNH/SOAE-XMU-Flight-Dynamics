function [ru] = My_density(H)

    if H < 11000
        temp = 1 - 0.0225569 * H / 1000;
        ru = 0.12492 * 9.8 * exp(4.255277 * log(temp));
    else
        temp = 11 - H / 1000;
        ru = 0.03718 * 9.8 * exp(temp / 6.318);
    end

end
