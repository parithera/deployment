#!/bin/bash

# Verify the presence of required commands
if ! command -v curl &>/dev/null; then
    echo "Error: curl is not installed."
    exit 1
fi

if ! command -v openssl &>/dev/null; then
    echo "Error: openssl is not installed."
    exit 1
fi

if ! command -v docker &>/dev/null || ! command -v docker compose &>/dev/null; then
    echo "Error: docker or docker compose is not installed."
    exit 1
fi

echo "All required commands are installed."

if [ -f ".git" ] || [ -d ".git" ]; then
    echo ".git found. Skipping repository clone."
else
    git clone https://github.com/parithera/deployment.git
    cd deployment
fi

# openssl needs to be installed
if [ ! -f jwt/private.pem ] && [ ! -f jwt/public.pem ]; then
    echo "Generating certificates"
    mkdir -p jwt
    openssl ecparam -name secp521r1 -genkey -noout -out jwt/private.pem
    openssl ec -in jwt/private.pem -pubout -out jwt/public.pem
else
    echo "Certificates found. Skipping certificate generation."
fi

read -p "Is this installation running on localhost (Y/n)? " local_install

domain_name="localhost"
if [[ "$local_install" == "n" || "$local_install" == "N" ]]; then
    read -p "Enter the domain name (localtest.io): " domain_name
    echo "Running installation for domain: $domain_name"
    sed -i.bak "s/WEB_HOST=https:\/\/localhost/WEB_HOST=https:\/\/$domain_name/" .env.api
    sed -i.bak "s/SERVER_NAME=localhost/SERVER_NAME=$domain_name/" .env.frontend
    read -p "Do you want Caddy to generate certificates (Y/n)? " caddy_generate_certs

    if [[ "$caddy_generate_certs" == "n" || "$caddy_generate_certs" == "N" ]]; then
        echo "Generating self-signed certificates..."
        mkdir -p certs
        openssl req -x509 -newkey rsa:4096 -keyout certs/tls.key -out certs/tls.pem -days 365 -nodes -subj "/CN=$domain_name"
        sed -i.bak 's/# - .\/certs:\/etc\/caddy\/certs:ro/- .\/certs:\/etc\/caddy\/certs:ro/' docker-compose.yaml
        sed -i.bak 's/# CADDY_SERVER_EXTRA_DIRECTIVES=/CADDY_SERVER_EXTRA_DIRECTIVES=/' .env.frontend
    else
        echo "Letting Caddy generate certificates..."
    fi
else
    echo "Running local installation."
fi

# Start DB container
docker compose -f docker-compose.yaml up db -d
# Download dumps
sleep 10
# Create Postgre databases
docker compose -f docker-compose.yaml -f docker-compose.knowledge.yaml run --rm knowledge -knowledge -action setup
# Restore database content from dumps
for db in "codeclarity" "knowledge" "config"; do
    docker compose -f docker-compose.yaml exec -T db sh -c "pg_restore -l /dump/$db.dump > /dump/$db.list && pg_restore -U postgres -d $db -L /dump/$db.list /dump/$db.dump"
done
# Start all containers
docker compose -f docker-compose.yaml up -d

echo "Installation successful, you can now visit: https://$domain_name:443"
